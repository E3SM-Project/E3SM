#!/usr/bin/env perl

use strict;
use warnings;
use XML::LibXML;

if ($ARGV[0] =~ '^(-)?-h(elp)?$') {
    die <<EOF;
SYNOPSIS
    compilers_translation_tool.pl INPUT_FILE OUTPUT_FILE
DESCRIPTION
    Converts an old-style config_compilers.xml file (INPUT_FILE) to a new-style
    config_build.xml file (OUTPUT_FILE). The new file should be checked both for
    the accuracy of comments, and the accuracy of <env>/<var> tags.

    The output file should conform to the config_build schema. This can be
    confirmed by running:

    xmllint -noout -schema \
        \$CIME_ROOT/config/xml_schemas/config_build.xsd OUTPUT_FILE
EOF
}

my $orig_doc = XML::LibXML->load_xml(location => $ARGV[0]);

my $orig_root = $orig_doc->documentElement();

# ==============================================================================
#
# Utility code
#
# ==============================================================================

# Returns false if both elements have an attribute with the given name, and
# the value of that attribute is not the same in the two cases. Otherwise
# returns true.
sub compatible_attribute {
    my ($element1, $element2, $attribute_name) = @_;
    if ($element1->hasAttribute($attribute_name) and
        $element2->hasAttribute($attribute_name)) {
        if ($element1->getAttribute($attribute_name) ne
            $element2->getAttribute($attribute_name)) {
            return 0;
        }
    }
    return 1;
}

# Returns true only if the MACH, OS, and COMPILER attributes are consistent
# between two elements.
sub machine_compatible_elements {
    my ($element1, $element2) = @_;
    return (compatible_attribute($element1, $element2, "MACH") and
            compatible_attribute($element1, $element2, "OS") and
            compatible_attribute($element1, $element2, "COMPILER"));
}

# Split up a text node by changing it into a text node, followed by an element,
# followed by another text node. This is used principly to replace non-XML
# syntax with XML-compatible markup.
#
# The return value is the three nodes that were created/modified within the
# parent of the input node.
sub split_text_node {
    my ($text_node, $before_text, $elem_name, $elem_text, $after_text) = @_;

    # What's the parent that we're replacing text in?
    my $parent_node = $text_node->parentNode;

    # In the child node, just keep the initial text.
    $text_node->setData($before_text);

    # Then, add a new element.
    my $element_text_node = XML::LibXML::Text->new($elem_text);
    my $element = XML::LibXML::Element->new($elem_name);
    $element->appendChild($element_text_node);
    $parent_node->insertAfter($element, $text_node);

    # Finally, add the remaining text.
    my $after_text_node = XML::LibXML::Text->new($after_text);
    $parent_node->insertAfter($after_text_node, $element);

    return ($text_node, $element_text_node, $after_text_node);
}

# ==============================================================================
#
# PHASE 1: Format tweaks
#
# Here we change a few things to conform to the new XML spec, but without doing
# any major restructuring of the XML tree.
#
# ==============================================================================

# MPILIB doesn't really need to exist on the "compiler" level, so push this
# attribute down to child nodes.
my @mpilib_compiler_nodes = $orig_root->findnodes("compiler[\@MPILIB]");

foreach my $mpilib_node ( @mpilib_compiler_nodes ) {
    # Get the value of the MPILIB attribute and remove it from the parent.
    my $mpilib = $mpilib_node->getAttribute("MPILIB");
    $mpilib_node->removeAttribute("MPILIB");

    # Add the MPILIB attribute to the children.
    my @child_nodes = $mpilib_node->childNodes();
    foreach my $child_node ( @child_nodes ) {
        # Only do this for element nodes.
        if (XML_ELEMENT_NODE == $child_node->nodeType) {
            # If the node already has an MPILIB attribute, then either the
            # attribute is the same as what we're going to add (and thus we
            # don't need to act), or it is different (and thus the child is
            # non-functional and we should remove it).
            if ($child_node->hasAttribute("MPILIB")) {
                my $child_mpilib = $child_node->getAttribute("MPILIB");
                if ($child_mpilib ne $mpilib) {
                    print STDERR "WARNING: Removing an element that had MPILIB=$child_mpilib\n";
                    print STDERR "This would not be used because the parent had MPILIB=$mpilib\n";
                    $mpilib_node->removeChild($child_node);
                }
            } else {
                $child_node->setAttribute("MPILIB", $mpilib);
            }
        }
    }
}

# Now we need to convert all the variable references in the file (these are
# mostly in the form "$(VAR_NAME)"), as well as shell commands.

# First, define a function to handle text nodes that contain variable
# references. Note that the first argument is the text node, while the
# second argument is the element immediately contained within the
# "compiler" element (which may not be the immediate parent for variable
# references within a shell command).
sub substitute_variable_references {
    my ($text_node, $variable_node) = @_;

    # Here we want to translate text of the form "$(FOO)" or "${FOO}" into
    # "<env>FOO</env>" or "<var>FOO</var>". But we have a problem, which is that
    # we can't tell whether a variable is supposed to be defined in a
    # config_compilers entry or in the environment. Therefore, we guess using a
    # few rules:
    #
    #  1) If there is no "FOO" defined in config_compilers, for the same
    #     compiler/OS/machine settings as the reference then we can be sure that
    #     it is an environment variable being referenced.
    #  2) If the reference is trivially self-referential (e.g. TRILINOS_PATH is
    #     set from TRILINOS_PATH), it must also be referring to an environment
    #     variable.
    #  3) If we have anything of the form "$ENV{FOO}", this is definitely an
    #     environment variable used by CMake.
    #  4) Otherwise, we assume that the variable is supposed to get its value
    #     from the config_compilers variable of the same name.
    my $before_text;
    my $variable_name;
    my $after_text;
    my $force_env = 0;
    my $variable_text = $text_node->getData();
    if ($variable_text =~ "(.*)\\\$\\\(([^)]*)\\\)(.*)") {
        $before_text = $1;
        $variable_name = $2;
        $after_text = $3;
    } elsif ($variable_text =~ "(.*)\\\$\\\{([^)]*)\\\}(.*)") {
        $before_text = $1;
        $variable_name = $2;
        $after_text = $3;
    } elsif ($variable_text =~ "(.*)\\\$ENV\\\{([^)]*)\\\}(.*)") {
        $before_text = $1;
        $variable_name = $2;
        $after_text = $3;
        $force_env = 1;
    } else {
        return;
    }

    # Default to assuming "env".
    my $element_name = "env";

    # Check if there's a variable in config_compilers itself that we
    # could be referring to, defined for the same machine.
    my $var_parent = $variable_node->parentNode;
    my @referenced_nodes = $orig_root->findnodes("compiler/$variable_name");
    foreach my $referenced_node (@referenced_nodes) {
        my $ref_parent = $referenced_node->parentNode;
        if (machine_compatible_elements($ref_parent, $var_parent)) {
            $element_name = "var";
            last;
        }
    }

    if (($variable_name eq $variable_node->nodeName) or $force_env) {
        $element_name = "env";
    }

    my ($before_node, $element_text_node, $after_node) =
        split_text_node($text_node, $before_text, $element_name,
                        $variable_name, $after_text);

    # There may be other variable references in the text, so now get those. We
    # can skip the variable name, though, since we don't allow nested variable
    # names in the above regex.
    substitute_variable_references($before_node, $variable_node);
    substitute_variable_references($after_node, $variable_node);
}


# Get all of the variables so that we can look at their contents.
my @variable_nodes = $orig_root->findnodes("compiler/*");

# Start with shell commands.
foreach my $variable_node ( @variable_nodes ) {
    my @child_nodes = $variable_node->childNodes();
    foreach my $child_node ( @child_nodes ) {
        # Only want to do something for (non-comment) text nodes.
        if (XML_TEXT_NODE == $child_node->nodeType) {
            my $variable_text = $child_node->getData();
            # This is a super-ugly regex, so I'll break it down:
            #  1. Capture everything before the interesting part.
            #  2. Find the string "$(shell ", then capture from here on.
            #  3. It can be followed by any number of characters that are not
            #     parentheses.
            #  4. There may also be a few sets of balanced parentheses, e.g.
            #     due to using things like "$(NETCDF_PATH)" in a shell command.
            #     We assume that there's no nesting, because in fact it is
            #     impossible to match against arbitrarily nested balanced
            #     parentheses without some amazing perl-specific regex-foo.
            #  5. There may be some number of characters that are not
            #     parentheses at the end of the shell command.
            #  6. Stop capturing the command right before the last ")".
            #  7. Capture everything after the interesting part.
            #
            #                      1   2               3     4                5       6    7
            if ($variable_text =~ "(.*)\\\$\\\(shell\ ([^()]*(?:\\\([^()]*\\\)[^()]*)?)\\\)(.*)") {
                my $before_text = $1;
                my $shell_command = $2;
                my $after_text = $3;

                my ($before_node, $element_text_node, $after_node) =
                    split_text_node($child_node, $before_text, "shell",
                                    $shell_command, $after_text);
                # There may be other variable references in the text, so now get those.
                substitute_variable_references($before_node, $variable_node);
                substitute_variable_references($element_text_node, $variable_node);
                substitute_variable_references($after_node, $variable_node);
            } else {
                substitute_variable_references($child_node, $variable_node);
            }
        }
    }
}

# Get rid of GPTL_CPPDEFS; we don't really to separate these from the regular
# CPPDEFS if we use MODEL="gptl".
my @gptl_cppdef_nodes = $orig_root->findnodes("compiler/GPTL_CPPDEFS");

foreach my $gptl_node ( @gptl_cppdef_nodes ) {
    # Even if this was just "GPTL_CPPDEFS" before, we want an "ADD_" prefix
    # because we don't want to overwrite other CPPDEFS.
    $gptl_node->setNodeName("ADD_CPPDEFS");
    $gptl_node->setAttribute("MODEL", "gptl");
}

# Same for the "ADD_" version.
@gptl_cppdef_nodes = $orig_root->findnodes("compiler/ADD_GPTL_CPPDEFS");

foreach my $gptl_node ( @gptl_cppdef_nodes ) {
    $gptl_node->setNodeName("ADD_CPPDEFS");
    $gptl_node->setAttribute("MODEL", "gptl");
}


# ==============================================================================
#
# PHASE 2: Revision
#
# This is the bulk of the work necessary to convert the document to conform to
# the new schema. We create a new document object, and load it up with
# information taken from the old one.
#
# ==============================================================================

# New XML document and root element.
my $doc = XML::LibXML::Document->new("1.0", "UTF-8");
my $root = $doc->createElement("config_build");
$doc->setDocumentElement($root);

# Before proceeding farther, there are a bunch of functions that we want to
# define.

# Just as a cleanup measure, let's sort elements as we go through the document.
# So we want a comparison operator to tell us which element goes first.
# First define a subroutine to compare on a particular attribute.
sub compare_element_on_attribute {
    my ($element1, $element2, $attribute) = @_;

    # Rules for this are:
    #  1. If both have the attribute, we rank them on that attribute with cmp.
    #  2. If one element has the attribute and the other doesn't, the one
    #     without the attribute comes first.
    #  3. If both lack the attribute, they are "equal" for our purposes.
    my $attr1 = $element1->getAttribute($attribute);
    my $attr2 = $element2->getAttribute($attribute);
    if (defined $attr1) {
        if (defined $attr2) {
            return $attr1 cmp $attr2;
        } else {
            return 1;
        }
    } else {
        if (defined $attr2) {
            return -1;
        } else {
            return 0;
        }
    }
}

# Now define an operator to do overall comparison of compiler elements.
sub compare_compiler_elements {
    my ($element1, $element2) = @_;

    # Sort by machine first.
    my $comp = compare_element_on_attribute($element1, $element2, "MACH");
    if ($comp != 0) {
        return $comp;
    }

    # Then sort by OS (this is mostly to ensure that OS-only settings come
    # before machine-specific settings overall).
    $comp = compare_element_on_attribute($element1, $element2, "OS");
    if ($comp != 0) {
        return $comp;
    }

    # Within each machine/OS, sort by compiler.
    return compare_element_on_attribute($element1, $element2, "COMPILER");
}

# This is for inserting whitespace to make the document more human-readable.
sub newline_node {
    return XML::LibXML::Text->new("\n");
}
sub two_spaces_node {
    return XML::LibXML::Text->new("  ");
}

# List of compiler flag elements (or items with similar format).
my %flag_element = map {$_, 1} ("CFLAGS", "CMAKE_OPTS", "CONFIG_ARGS",
                     "CPPDEFS", "CXX_LDFLAGS", "CXX_LIBS", "FC_AUTO_R8",
                     "FFLAGS", "FFLAGS_NOOPT", "FIXEDFLAGS", "FREEFLAGS",
                     "LDFLAGS", "MLIBS", "SLIBS");

# All other elements.
my %non_flag_element = map {$_, 1} ("ALBANY_PATH", "CONFIG_SHELL", "CPRE",
                         "CXX_LINKER", "ESMF_LIBDIR", "HAS_F2008_CONTIGUOUS",
                         "HDF5_PATH", "LAPACK_LIBDIR", "LD", "MPI_LIB_NAME",
                         "MPI_PATH", "MPICC", "MPICXX", "MPIFC", "NETCDF_PATH",
                         "PAPI_INC", "PAPI_LIB", "PETSC_PATH", "PFUNIT_PATH",
                         "PIO_FILESYSTEM_HINTS", "PNETCDF_PATH", "SCC", "SCXX",
                         "SFC", "SUPPORTS_CXX", "TRILINOS_PATH");

# Function to merge a compiler element from the old format into the new one.
sub merge_element {
    my ($orig_node, $new_node) = @_;

    my @comment_stack;
    foreach my $orig_child ($orig_node->childNodes()) {
        if ($orig_child->nodeType == XML_ELEMENT_NODE) {
            # Name, and basename with "ADD_" stripped off.
            my $name = $orig_child->nodeName;
            my $basename = $name =~ s/^ADD_//r;
            if ($non_flag_element{$name}) {
                # Use "1" to get a deep clone, copying all child nodes.
                my $new_child = $orig_child->cloneNode(1);
                # Find out where in the new document we should put the element.
                # This is just using alphabetical order.
                my $prev_node = $new_node->firstChild;
                foreach my $place_node ($new_node->childNodes()) {
                    unless ($place_node->nodeType == XML_ELEMENT_NODE) {next};
                    if ($place_node->nodeName gt $name) {
                        # If we've reached a node that should be "below" this
                        # one, put in our new child node here.
                        last;
                    }
                    # This should get the newline after this node.
                    $prev_node = $place_node->nextSibling();
                }
                $new_node->insertAfter(newline_node(), $prev_node);
                $new_node->insertAfter($new_child, $prev_node);
                $new_node->insertAfter(two_spaces_node(), $prev_node);
                # Deal with comments internal to the block.
                my $comment = pop @comment_stack;
                while (defined $comment) {
                    $new_node->insertAfter(newline_node(), $prev_node);
                    $new_node->insertAfter($comment, $prev_node);
                    $new_node->insertAfter(two_spaces_node(), $prev_node);
                    $comment = pop @comment_stack;
                }
            } elsif ($flag_element{$basename}) {
                # First ensure that a block for this set of flags exists at all.
                my @flag_blocks = $new_node->findnodes($basename);
                if ((scalar @flag_blocks) > 1) {
                    die "Somehow have more than one $basename element for the same compiler.";
                }
                # Get the block if it exists.
                my $flag_block = $flag_blocks[0];
                # If not, do something similar to the above (slightly simplified
                # because we aren't dealing with comments.
                if (! defined $flag_block) {
                    # Create a new element.
                    my $new_child = XML::LibXML::Element->new($basename);
                    $new_child->appendChild(newline_node());
                    $new_child->appendChild(two_spaces_node());
                    # Find out where in the new document we should put the element.
                    # This is just using alphabetical order.
                    my $prev_node = $new_node->firstChild;
                    foreach my $place_node ($new_node->childNodes()) {
                        unless ($place_node->nodeType == XML_ELEMENT_NODE) {next};
                        if ($place_node->nodeName gt $basename) {
                            # If we've reached a node that should be "below" this
                            # one, put in our new child node here.
                            last;
                        }
                        # This should get the newline after this node.
                        $prev_node = $place_node->nextSibling();
                    }
                    $new_node->insertAfter(newline_node(), $prev_node);
                    $new_node->insertAfter($new_child, $prev_node);
                    $new_node->insertAfter(two_spaces_node(), $prev_node);
                    $flag_block = $new_child;
                }
                # Now that we have the block, add the information to it.
                # Start by cloning the original node.
                my $new_child = $orig_child->cloneNode(1);
                # If the original node name had an "ADD_" prefix, we want to
                # create an "append" tag, otherwise a "base" tag.
                if ($name eq $basename) {
                    $new_child->setNodeName("base");
                    $flag_block->insertAfter(newline_node(), $flag_block->firstChild);
                    $flag_block->insertAfter($new_child, $flag_block->firstChild);
                    $flag_block->insertAfter(two_spaces_node(), $flag_block->firstChild);
                    $flag_block->insertAfter(two_spaces_node(), $flag_block->firstChild);
                    my $comment = pop @comment_stack;
                    while (defined $comment) {
                        $flag_block->insertAfter(newline_node(), $flag_block->firstChild);
                        $flag_block->insertAfter($comment, $flag_block->firstChild);
                        $flag_block->insertAfter(two_spaces_node(), $flag_block->firstChild);
                        $flag_block->insertAfter(two_spaces_node(), $flag_block->firstChild);
                        $comment = pop @comment_stack;
                    }
                } else {
                    my $comment = shift @comment_stack;
                    while (defined $comment) {
                        $flag_block->appendChild(two_spaces_node());
                        $flag_block->appendChild($comment);
                        $flag_block->appendChild(newline_node());
                        $flag_block->appendChild(two_spaces_node());
                        $comment = shift @comment_stack;
                    }
                    $new_child->setNodeName("append");
                    $flag_block->appendChild(two_spaces_node());
                    $flag_block->appendChild($new_child);
                    $flag_block->appendChild(newline_node());
                    $flag_block->appendChild(two_spaces_node());
                }
            } else {
                die "Unrecognized variable in compiler element: $name";
            }
        } elsif ($orig_child->nodeType == XML_COMMENT_NODE) {
            # Comment nodes should be copied to the new file. However, since
            # they usually refer to whatever comes next, we want to make a stack
            # that we will push into the new document later.
            push @comment_stack, XML::LibXML::Comment->new($orig_child->data);
        }
        # Anything that isn't a comment or element node is probably whitespace
        # and can be ignored.
    }
}

# Function to actually create a new compiler element.
sub new_format_element {
    my ($orig_node) = @_;

    # Make a new element, and add appropriate attributes from the original node.
    my $new_node = XML::LibXML::Element->new("compiler");
    foreach my $attribute ("MACH", "OS", "COMPILER") {
        my $attr_value = $orig_node->getAttribute($attribute);
        # If the original element has an attribute, add it to
        # the new element.
        if (defined $attr_value) {
            $new_node->setAttribute($attribute, $attr_value);
        }
    }

    # Start with a newline.
    $new_node->appendChild(newline_node());

    merge_element($orig_node, $new_node);

    return $new_node;
}

# Now that we have all these functions, proceed to create the new document.

# Add a new line to the document right up front.
$root->appendChild(newline_node());

my @comment_stack;

# Migrating all compiler and comment nodes from the original document.
foreach my $orig_node ($orig_root->childNodes()) {
    # If we have an element node, we need to handle "compiler" elements.
    if ($orig_node->nodeType == XML_ELEMENT_NODE) {
        if ($orig_node->nodeName eq "compiler") {
            # Find out where in the new document we should put the element.
            my $element_placed = 0;
            foreach my $place_node ($root->findnodes("compiler")) {
                my $elem_comp = compare_compiler_elements($orig_node, $place_node);
                if ($elem_comp == -1) {
                    # If the element we're adding comes before this item, place
                    # the new element now.

                    # Deal with comment nodes first.
                    my $comment = shift @comment_stack;
                    while (defined $comment) {
                        $root->insertBefore($comment, $place_node);
                        $root->insertBefore(newline_node(), $place_node);
                        $comment = shift @comment_stack;
                    }
                    # Create compiler element.
                    my $element = new_format_element($orig_node);

                    # Place it, and put in some whitespace.
                    $root->insertBefore($element, $place_node);
                    $root->insertBefore(newline_node(), $place_node);
                    $root->insertBefore(newline_node(), $place_node);
                    $element_placed = 1;
                    last;
                } elsif ($elem_comp == 0) {
                    # If a matching element is present, don't add another one.
                    # Just merge the element in place and quit.
                    merge_element($orig_node, $place_node);
                    $element_placed = 1;
                    last;
                }
                # If the element we're adding comes after this item, just keep
                # going.
            }
            # If we didn't place the element, do so at the end of the list.
            if ($element_placed == 0) {
                # Deal with comment nodes first.
                my $comment = shift @comment_stack;
                while (defined $comment) {
                    $root->appendChild($comment);
                    $root->appendChild(newline_node());
                    $comment = shift @comment_stack;
                }
                # Create compiler element.
                my $element = new_format_element($orig_node);

                # Place it, and put in some whitespace.
                $root->appendChild($element);
                $root->appendChild(newline_node());
                $root->appendChild(newline_node());
            }
        } else {
            my $bad_name = $orig_node->nodeName;
            die "Unrecognized node at root level: $bad_name";
        }
    } elsif ($orig_node->nodeType == XML_COMMENT_NODE) {
        # Comment nodes should be copied to the new file. However, since they
        # usually refer to whatever comes next, we want to make a stack that we
        # will push into the new document later.
        push @comment_stack, XML::LibXML::Comment->new($orig_node->data);
    }
    # Anything that isn't a comment or element node is probably whitespace
    # and can be ignored.
}

# Print the completed product!
$doc->toFile($ARGV[1]);

print STDERR "WARNING: It is not always possible to tell whether certain \n".
    "variables, e.g. NETCDF_PATH, are intended to come from the environment \n".
    "or from config_build itself. Please check all <var> tags to ensure \n".
    "that they should not be changed to <env> instead.\n";

print STDERR "WARNING: Comments have been copied over from the old file. \n".
    "However, there is no way to be certain that the comments are still in \n".
    "the right place in the XML file, nor to be sure that they even still \n".
    "make sense. The user is encouraged to review all comments in the \n".
    "generated file for accuracy.\n";
