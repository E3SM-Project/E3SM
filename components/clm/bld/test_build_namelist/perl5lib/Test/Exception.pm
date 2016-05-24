use strict;
use warnings;

package Test::Exception;
use Test::Builder;
use Sub::Uplevel qw( uplevel );
use base qw( Exporter );

our $VERSION = '0.32';
our @EXPORT = qw(dies_ok lives_ok throws_ok lives_and);

my $Tester = Test::Builder->new;

sub import {
    my $self = shift;
    if ( @_ ) {
        my $package = caller;
        $Tester->exported_to( $package );
        $Tester->plan( @_ );
    };
    $self->export_to_level( 1, $self, $_ ) foreach @EXPORT;
}

=head1 NAME

Test::Exception - Test exception based code

=head1 SYNOPSIS

  use Test::More tests => 5;
  use Test::Exception;

  # or if you don't need Test::More

  use Test::Exception tests => 5;

  # then...

  # Check that the stringified exception matches given regex
  throws_ok { $foo->method } qr/division by zero/, 'zero caught okay';

  # Check an exception of the given class (or subclass) is thrown
  throws_ok { $foo->method } 'Error::Simple', 'simple error thrown';
  
  # all Test::Exceptions subroutines are guaranteed to preserve the state 
  # of $@ so you can do things like this after throws_ok and dies_ok
  like $@, 'what the stringified exception should look like';

  # Check that something died - we do not care why
  dies_ok { $foo->method } 'expecting to die';

  # Check that something did not die
  lives_ok { $foo->method } 'expecting to live';

  # Check that a test runs without an exception
  lives_and { is $foo->method, 42 } 'method is 42';
  
  # or if you don't like prototyped functions
  
  throws_ok( sub { $foo->method }, qr/division by zero/,
      'zero caught okay' );
  throws_ok( sub { $foo->method }, 'Error::Simple', 
      'simple error thrown' );
  dies_ok( sub { $foo->method }, 'expecting to die' );
  lives_ok( sub { $foo->method }, 'expecting to live' );
  lives_and( sub { is $foo->method, 42 }, 'method is 42' );


=head1 DESCRIPTION

This module provides a few convenience methods for testing exception based code. It is built with 
L<Test::Builder> and plays happily with L<Test::More> and friends.

If you are not already familiar with L<Test::More> now would be the time to go take a look.

You can specify the test plan when you C<use Test::Exception> in the same way as C<use Test::More>.
See L<Test::More> for details.

NOTE: Test::Exception only checks for exceptions. It will ignore other methods of stopping 
program execution - including exit(). If you have an exit() in evalled code Test::Exception
will not catch this with any of its testing functions.

=cut

sub _quiet_caller (;$) { ## no critic Prototypes
    my $height = $_[0];
    $height++;

    if ( CORE::caller() eq 'DB' ) {
        # passthrough the @DB::args trick
        package DB;
        if( wantarray ) {
            if ( !@_ ) {
                return (CORE::caller($height))[0..2];
            }
            else {
                # If we got here, we are within a Test::Exception test, and
                # something is producing a stacktrace. In case this is a full
                # trace (i.e. confess() ), we have to make sure that the sub
                # args are not visible. If we do not do this, and the test in
                # question is throws_ok() with a regex, it will end up matching
                # against itself in the args to throws_ok().
                #
                # While it is possible (and maybe wise), to test if we are
                # indeed running under throws_ok (by crawling the stack right
                # up from here), the old behavior of Test::Exception was to
                # simply obliterate @DB::args altogether in _quiet_caller, so
                # we are just preserving the behavior to avoid surprises
                #
                my @frame_info = CORE::caller($height);
                @DB::args = ();
                return @frame_info;
            }
        }

        # fallback if nothing above returns
        return CORE::caller($height);
    }
    else {
        if( wantarray and !@_ ) {
            return (CORE::caller($height))[0..2];
        }
        else {
            return CORE::caller($height);
        }
    }
}

sub _try_as_caller {
    my $coderef = shift;

    # local works here because Sub::Uplevel has already overridden caller
    local *CORE::GLOBAL::caller;
    { no warnings 'redefine'; *CORE::GLOBAL::caller = \&_quiet_caller; }

    eval { uplevel 3, $coderef };
    return $@;
};


sub _is_exception {
    my $exception = shift;
    return ref $exception || $exception ne '';
};


sub _exception_as_string {
    my ( $prefix, $exception ) = @_;
    return "$prefix normal exit" unless _is_exception( $exception );
    my $class = ref $exception;
    $exception = "$class ($exception)" 
            if $class && "$exception" !~ m/^\Q$class/;
    chomp $exception;
    return "$prefix $exception";
};


=over 4

=item B<throws_ok>

Tests to see that a specific exception is thrown. throws_ok() has two forms: 

  throws_ok BLOCK REGEX, TEST_DESCRIPTION
  throws_ok BLOCK CLASS, TEST_DESCRIPTION

In the first form the test passes if the stringified exception matches the give regular expression. For example:

    throws_ok { read_file( 'unreadable' ) } qr/No file/, 'no file';

If your perl does not support C<qr//> you can also pass a regex-like string, for example:

    throws_ok { read_file( 'unreadable' ) } '/No file/', 'no file';

The second form of throws_ok() test passes if the exception is of the same class as the one supplied, or a subclass of that class. For example:

    throws_ok { $foo->bar } "Error::Simple", 'simple error';

Will only pass if the C<bar> method throws an Error::Simple exception, or a subclass of an Error::Simple exception.

You can get the same effect by passing an instance of the exception you want to look for. The following is equivalent to the previous example:

    my $SIMPLE = Error::Simple->new;
    throws_ok { $foo->bar } $SIMPLE, 'simple error';

Should a throws_ok() test fail it produces appropriate diagnostic messages. For example:

    not ok 3 - simple error
    #     Failed test (test.t at line 48)
    # expecting: Error::Simple exception
    # found: normal exit

Like all other Test::Exception functions you can avoid prototypes by passing a subroutine explicitly:

    throws_ok( sub {$foo->bar}, "Error::Simple", 'simple error' );

A true value is returned if the test succeeds, false otherwise. On exit $@ is guaranteed to be the cause of death (if any).

A description of the exception being checked is used if no optional test description is passed.

NOTE: Rememeber when you C<die $string_without_a_trailing_newline> perl will 
automatically add the current script line number, input line number and a newline. This will
form part of the string that throws_ok regular expressions match against.


=cut


sub throws_ok (&$;$) {
    my ( $coderef, $expecting, $description ) = @_;
    unless (defined $expecting) {
      require Carp;
      Carp::croak( "throws_ok: must pass exception class/object or regex" ); 
    }
    $description = _exception_as_string( "threw", $expecting )
    	unless defined $description;
    my $exception = _try_as_caller( $coderef );
    my $regex = $Tester->maybe_regex( $expecting );
    my $ok = $regex 
        ? ( $exception =~ m/$regex/ ) 
        : eval { 
            $exception->isa( ref $expecting ? ref $expecting : $expecting ) 
        };
    $Tester->ok( $ok, $description );
    unless ( $ok ) {
        $Tester->diag( _exception_as_string( "expecting:", $expecting ) );
        $Tester->diag( _exception_as_string( "found:", $exception ) );
    };
    $@ = $exception;
    return $ok;
};


=item B<dies_ok>

Checks that a piece of code dies, rather than returning normally. For example:

    sub div {
        my ( $a, $b ) = @_;
        return $a / $b;
    };

    dies_ok { div( 1, 0 ) } 'divide by zero detected';

    # or if you don't like prototypes
    dies_ok( sub { div( 1, 0 ) }, 'divide by zero detected' );

A true value is returned if the test succeeds, false otherwise. On exit $@ is guaranteed to be the cause of death (if any).

Remember: This test will pass if the code dies for any reason. If you care about the reason it might be more sensible to write a more specific test using throws_ok().

The test description is optional, but recommended. 

=cut

sub dies_ok (&;$) {
    my ( $coderef, $description ) = @_;
    my $exception = _try_as_caller( $coderef );
    my $ok = $Tester->ok( _is_exception($exception), $description );
    $@ = $exception;
    return $ok;
}


=item B<lives_ok>

Checks that a piece of code doesn't die. This allows your test script to continue, rather than aborting if you get an unexpected exception. For example:

    sub read_file {
        my $file = shift;
        local $/;
        open my $fh, '<', $file or die "open failed ($!)\n";
        $file = <FILE>;
        return $file;
    };

    my $file;
    lives_ok { $file = read_file('test.txt') } 'file read';

    # or if you don't like prototypes
    lives_ok( sub { $file = read_file('test.txt') }, 'file read' );

Should a lives_ok() test fail it produces appropriate diagnostic messages. For example:

    not ok 1 - file read
    #     Failed test (test.t at line 15)
    # died: open failed (No such file or directory)

A true value is returned if the test succeeds, false otherwise. On exit $@ is guaranteed to be the cause of death (if any).

The test description is optional, but recommended. 

=cut

sub lives_ok (&;$) {
    my ( $coderef, $description ) = @_;
    my $exception = _try_as_caller( $coderef );
    my $ok = $Tester->ok( ! _is_exception( $exception ), $description );
	$Tester->diag( _exception_as_string( "died:", $exception ) ) unless $ok;
    $@ = $exception;
    return $ok;
}


=item B<lives_and>

Run a test that may throw an exception. For example, instead of doing:

  my $file;
  lives_ok { $file = read_file('answer.txt') } 'read_file worked';
  is $file, "42", 'answer was 42';

You can use lives_and() like this:

  lives_and { is read_file('answer.txt'), "42" } 'answer is 42';
  # or if you don't like prototypes
  lives_and(sub {is read_file('answer.txt'), "42"}, 'answer is 42');

Which is the same as doing

  is read_file('answer.txt'), "42\n", 'answer is 42';

unless C<read_file('answer.txt')> dies, in which case you get the same kind of error as lives_ok()

  not ok 1 - answer is 42
  #     Failed test (test.t at line 15)
  # died: open failed (No such file or directory)

A true value is returned if the test succeeds, false otherwise. On exit $@ is guaranteed to be the cause of death (if any).

The test description is optional, but recommended.

=cut

sub lives_and (&;$) {
    my ( $test, $description ) = @_;
    {
        local $Test::Builder::Level = $Test::Builder::Level + 1;
        my $ok = \&Test::Builder::ok;
        no warnings;
        local *Test::Builder::ok = sub {
            $_[2] = $description unless defined $_[2];
            $ok->(@_);
        };
        use warnings;
        eval { $test->() } and return 1;
    };
    my $exception = $@;
    if ( _is_exception( $exception ) ) {
        $Tester->ok( 0, $description );
        $Tester->diag( _exception_as_string( "died:", $exception ) );
    };
    $@ = $exception;
    return;
}

=back


=head1 SKIPPING TEST::EXCEPTION TESTS

Sometimes we want to use Test::Exception tests in a test suite, but don't want to force the user to have Test::Exception installed. One way to do this is to skip the tests if Test::Exception is absent. You can do this with code something like this:

  use strict;
  use warnings;
  use Test::More;
  
  BEGIN {
      eval "use Test::Exception";
      plan skip_all => "Test::Exception needed" if $@;
  }
  
  plan tests => 2;
  # ... tests that need Test::Exception ...

Note that we load Test::Exception in a C<BEGIN> block ensuring that the subroutine prototypes are in place before the rest of the test script is compiled.


=head1 BUGS

There are some edge cases in Perl's exception handling where Test::Exception will miss exceptions
thrown in DESTROY blocks. See the RT bug L<http://rt.cpan.org/Ticket/Display.html?id=24678> for
details, along with the t/edge-cases.t in the distribution test suite. These will be addressed in
a future Test::Exception release.

If you find any more bugs please let me know by e-mail, or report the problem with 
L<http://rt.cpan.org/>.


=head1 COMMUNITY

=over 4

=item perl-qa

If you are interested in testing using Perl I recommend you visit L<http://qa.perl.org/> and join the excellent perl-qa mailing list. See L<http://lists.perl.org/showlist.cgi?name=perl-qa> for details on how to subscribe.

=item perlmonks

You can find users of Test::Exception, including the module author, on  L<http://www.perlmonks.org/>. Feel free to ask questions on Test::Exception there.

=item CPAN::Forum

The CPAN Forum is a web forum for discussing Perl's CPAN modules.   The Test::Exception forum can be found at L<http://www.cpanforum.com/dist/Test-Exception>.

=item AnnoCPAN

AnnoCPAN is a web site that allows community annotations of Perl module documentation. The Test::Exception annotations can be found at L<http://annocpan.org/~ADIE/Test-Exception/>.

=back


=head1 TO DO

If you think this module should do something that it doesn't (or does something that it shouldn't) please let me know.

You can see my current to do list at L<http://adrianh.tadalist.com/lists/public/15421>, with an RSS feed of changes at L<http://adrianh.tadalist.com/lists/feed_public/15421>.


=head1 ACKNOWLEDGMENTS

Thanks to chromatic and Michael G Schwern for the excellent Test::Builder, without which this module wouldn't be possible.

Thanks to 
Adam Kennedy,
Andy Lester, 
Aristotle Pagaltzis, 
Ben Prew, 
Cees Hek,
Chris Dolan,
chromatic, 
Curt Sampson,
David Cantrell,
David Golden, 
David Tulloh,
David Wheeler, 
J. K. O'Brien,
Janek Schleicher,
Jim Keenan, 
Jos I. Boumans, 
Joshua ben Jore,
Jost Krieger,
Mark Fowler, 
Michael G Schwern, 
Nadim Khemir,
Paul McCann,
Perrin Harkins, 
Peter Rabbitson,
Peter Scott, 
Ricardo Signes,
Rob Muhlestein,
Scott R. Godin,
Steve Purkis,
Steve, 
Tim Bunce,
and various anonymous folk for comments, suggestions, bug reports and patches.

=head1 AUTHOR

Adrian Howard <adrianh@quietstars.com>

If you can spare the time, please drop me a line if you find this module useful.


=head1 SEE ALSO

=over 4

=item L<http://del.icio.us/tag/Test::Exception>

Delicious links on Test::Exception.

=item L<Test::Warn> & L<Test::NoWarnings>

Modules to help test warnings.

=item L<Test::Builder>

Support module for building test libraries.

=item L<Test::Simple> & L<Test::More>

Basic utilities for writing tests.

=item L<http://qa.perl.org/test-modules.html>

Overview of some of the many testing modules available on CPAN.

=item L<http://del.icio.us/tag/perl+testing>

Delicious links on perl testing.

=back


=head1 LICENCE

Copyright 2002-2007 Adrian Howard, All Rights Reserved.

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut

1;
