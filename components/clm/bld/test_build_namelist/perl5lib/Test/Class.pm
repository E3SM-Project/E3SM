use strict;
use warnings;
use 5.006;

package Test::Class;

use Attribute::Handlers;
use Carp;
use MRO::Compat;
use Storable qw(dclone);
use Test::Builder;
use Test::Class::MethodInfo;

our $VERSION = '0.41';

my $Check_block_has_run;
{
    no warnings 'void';
    CHECK { $Check_block_has_run = 1 };
}

use constant NO_PLAN    => "no_plan";
use constant SETUP      => "setup";
use constant TEST       => "test";
use constant TEARDOWN   => "teardown";
use constant STARTUP    => "startup";
use constant SHUTDOWN   => "shutdown";


our $Current_method = undef;
sub current_method { $Current_method };


my $Builder = Test::Builder->new;
sub builder { $Builder };


my $Tests = {};
my @Filters = ();


my %_Test;  # inside-out object field indexed on $self

sub DESTROY {
    my $self = shift;
    delete $_Test{ $self };
};

sub _test_info {
    my $self = shift;
    return ref($self) ? $_Test{$self} : $Tests;
};

sub _method_info {
    my ($self, $class, $method) = @_;
    return( _test_info($self)->{$class}->{$method} );
};

sub _methods_of_class {
    my ( $self, $class ) = @_;
    my $test_info = _test_info($self) 
        or die "Test::Class internals seem confused. Did you override "
            . "new() in a sub-class or via multiple inheritence?\n";
    return values %{ $test_info->{$class} };
};

sub _parse_attribute_args {
    my $args = shift || '';
    my $num_tests;
    my $type;
    $args =~ s/\s+//sg;
    foreach my $arg (split /=>/, $args) {
        if (Test::Class::MethodInfo->is_num_tests($arg)) {
            $num_tests = $arg;
        } elsif (Test::Class::MethodInfo->is_method_type($arg)) {
            $type = $arg;
        } else {
            die 'bad attribute args';
        };
    };
    return( $type, $num_tests );
};

sub _is_public_method {
    my ($class, $name) = @_;
    my @parents = @{mro::get_linear_isa($class)};
    shift @parents;
    foreach my $parent_class ( @parents ) {
        return unless $parent_class->can( $name );
        return if _method_info( $class, $parent_class, $name );
    }
    return 1;
}

sub Test : ATTR(CODE,RAWDATA) {
    my ($class, $symbol, $code_ref, $attr, $args) = @_;
    if ($symbol eq "ANON") {
        warn "cannot test anonymous subs - you probably loaded a Test::Class too late (after the CHECK block was run). See 'A NOTE ON LOADING TEST CLASSES' in perldoc Test::Class for more details\n";
    } else {
        my $name = *{$symbol}{NAME};
        warn "overriding public method $name with a test method in $class\n"
                if _is_public_method( $class, $name );
        eval { $class->add_testinfo($name, _parse_attribute_args($args)) } 
            || warn "bad test definition '$args' in $class->$name\n";   
    };
};

sub Tests : ATTR(CODE,RAWDATA) {
    my ($class, $symbol, $code_ref, $attr, $args) = @_;
    $args ||= 'no_plan';
    Test( $class, $symbol, $code_ref, $attr, $args );
};

sub add_testinfo {
    my($class, $name, $type, $num_tests) = @_;
    $Tests->{$class}->{$name} = Test::Class::MethodInfo->new(
        name => $name,
        num_tests => $num_tests,
        type => $type,
    );
}

sub _class_of {
    my $self = shift;
    return ref $self ? ref $self : $self;
}

sub new {
    my $proto = shift;
    my $class = _class_of( $proto );
    $proto = {} unless ref($proto);
    my $self = bless {%$proto, @_}, $class;
    $_Test{$self} = dclone($Tests);
    return($self);
};

sub _get_methods {
    my ( $self, @types ) = @_;
    my $test_class = _class_of( $self );
    
    my $test_method_regexp = $ENV{ TEST_METHOD } || '.*';
    my $method_regexp = eval { qr/\A$test_method_regexp\z/ };
    die "TEST_METHOD ($test_method_regexp) is not a valid regexp: $@" if $@;
    
    my %methods = ();
    foreach my $class ( @{mro::get_linear_isa( $test_class )} ) {
      FILTER:
        foreach my $info ( _methods_of_class( $self, $class ) ) {
            my $name = $info->name;

            if ( $info->type eq TEST ) {
                # determine if method is filtered, true if *any* filter
                # returns false.
                foreach my $filter ( @Filters ) {
                    next FILTER unless $filter->( $class, $name );
                }
            }

            foreach my $type ( @types ) {
                if ( $info->is_type( $type ) ) {
                    $methods{ $name } = 1 
                        unless $type eq TEST && $name !~ $method_regexp;
                }
            };
        };
    };

    my @methods = sort keys %methods;
    return @methods;
};

sub _num_expected_tests {
    my $self = shift;
    if (my $reason = $self->SKIP_CLASS ) {
       return $reason eq "1" ? 0 : 1;
    };
    my @test_methods = _get_methods($self, TEST);
    return 0 unless @test_methods;
    my @startup_shutdown_methods = 
            _get_methods($self, STARTUP, SHUTDOWN);
    my $num_startup_shutdown_methods = 
            _total_num_tests($self, @startup_shutdown_methods);
    return(NO_PLAN) if $num_startup_shutdown_methods eq NO_PLAN;
    my @fixture_methods = _get_methods($self, SETUP, TEARDOWN);
    my $num_fixture_tests = _total_num_tests($self, @fixture_methods);
    return(NO_PLAN) if $num_fixture_tests eq NO_PLAN;
    my $num_tests = _total_num_tests($self, @test_methods);
    return(NO_PLAN) if $num_tests eq NO_PLAN;
    return($num_startup_shutdown_methods + $num_tests + @test_methods * $num_fixture_tests);
};

sub expected_tests {
    my $total = 0;
    foreach my $test (@_) {
        if ( _isa_class( __PACKAGE__, $test ) ) {
            my $n = _num_expected_tests($test);
            return NO_PLAN if $n eq NO_PLAN;
            $total += $n;
        } elsif ( defined $test && $test =~ m/^\d+$/ ) {
            $total += $test;
        } else {
            $test = 'undef' unless defined $test;
            croak "$test is not a Test::Class or an integer";
        };
    };
    return $total;
};

sub _total_num_tests {
    my ($self, @methods) = @_;
    my $class = _class_of( $self );
    my $total_num_tests = 0;
    foreach my $method (@methods) {
        foreach my $class (@{mro::get_linear_isa($class)}) {
            my $info = _method_info($self, $class, $method);
            next unless $info;
            my $num_tests = $info->num_tests;
            return(NO_PLAN) if ($num_tests eq NO_PLAN);
            $total_num_tests += $num_tests;
            last unless $num_tests =~ m/^\+/
        };
    };
    return($total_num_tests);
};

sub _has_no_tests {
    my ( $self, $method ) = @_;
    return _total_num_tests( $self, $method ) eq '0';
}

sub _all_ok_from {
    my ($self, $start_test) = @_;

    # The Test::Builder 1.5 way to do it
    if( $Builder->can("history") ) {
        return $Builder->history->can_succeed;
    }
    # The Test::Builder 0.x way to do it
    else {
        my $current_test = $Builder->current_test;
        return(1) if $start_test == $current_test;
        my @results = ($Builder->summary)[$start_test .. $current_test-1];
        foreach my $result (@results) { return(0) unless $result };
        return(1);
    }
};

sub _exception_failure {
    my ($self, $method, $exception, $tests) = @_;
    local $Test::Builder::Level = 3;
    my $message = $method;
    $message .= " (for test method '$Current_method')"
            if defined $Current_method && $method ne $Current_method;
    _show_header($self, @$tests);
    $Builder->ok(0, "$message died ($exception)");
    _threw_exception( $self, $method => 1 );
};

my %threw_exception;
sub _threw_exception {
    my ( $self, $method, $optional_value) = @_;
    my $class = ref( $self );
    $threw_exception{ $class }{ $method } = $optional_value 
        if defined $optional_value;
    return $threw_exception{ $class }{ $method };
}

sub _run_method {
    my ($self, $method, $tests) = @_;
    _threw_exception( $self, $method => 0 );
    my $num_start = $Builder->current_test;
    my $skip_reason;
    my $original_ok = \&Test::Builder::ok;
    no warnings;
    local *Test::Builder::ok = sub {
        my ($builder, $test, $description) = @_;
        local $Test::Builder::Level = $Test::Builder::Level+1;
        unless ( defined($description) ) {
            $description = $self->current_method;
            $description =~ tr/_/ /;
        };
        my $is_ok = $original_ok->($builder, $test, $description);
        unless ( $is_ok ) {
            my $class = ref $self;
            $Builder->diag( "  (in $class->$method)" );
        };
        return $is_ok;
    };
    $skip_reason = eval {$self->$method};
    $skip_reason = $method unless $skip_reason;
    my $exception = $@;
    chomp($exception) if $exception;
    my $num_done = $Builder->current_test - $num_start;
    my $num_expected = _total_num_tests($self, $method);
    $num_expected = $num_done if $num_expected eq NO_PLAN;
    if ($num_done == $num_expected) {
        _exception_failure($self, $method, $exception, $tests) 
                unless $exception eq '';
    } elsif ($num_done > $num_expected) {
        my $class = ref $self;
        $Builder->diag("expected $num_expected test(s) in $class\::$method, $num_done completed\n");
    } else {
        until (($Builder->current_test - $num_start) >= $num_expected) {
            if ($exception ne '') {
                _exception_failure($self, $method, $exception, $tests);
                $skip_reason = "$method died";
                $exception = '';
            } else {
                if ($self->fail_if_returned_early) {
                    my $class = ref $self;
                    $Builder->ok(0, "($class\::$method returned before plan complete)");
                } else {
                    $Builder->skip( $skip_reason );
                }
            };
        };
    };
    return(_all_ok_from($self, $num_start));
};

sub fail_if_returned_early { 0 }

sub _show_header {
    my ($self, @tests) = @_;
    return if $Builder->has_plan;
    my $num_tests = Test::Class->expected_tests(@tests);
    if ($num_tests eq NO_PLAN) {
        $Builder->no_plan;
    } else {
        $Builder->expected_tests($num_tests);
    };
};

my %SKIP_THIS_CLASS = ();

sub SKIP_CLASS {
    my $class = shift;
    $SKIP_THIS_CLASS{ $class } = shift if @_;
    return $SKIP_THIS_CLASS{ $class };
};

sub _isa_class {
    my ( $class, $object_or_class ) = @_;
    return unless defined $object_or_class;
    return if $object_or_class eq 'Contextual::Return::Value';
    return eval { 
        $object_or_class->isa( $class ) and $object_or_class->can( 'runtests' )
    };
}

sub _test_classes {
    my $class = shift;
    return( @{mro::get_isarev($class)}, $class );
};

sub runtests {
    die "Test::Class was loaded too late (after the CHECK block was run). See 'A NOTE ON LOADING TEST CLASSES' in perldoc Test::Class for more details\n"
        unless $Check_block_has_run;
    my @tests = @_;
    if (@tests == 1 && !ref($tests[0])) {
        my $base_class = shift @tests;
        @tests = _test_classes( $base_class );
    };
    my $all_passed = 1;
    TEST_OBJECT: foreach my $t (@tests) {
        # SHOULD ALSO ALLOW NO_PLAN
        next if $t =~ m/^\d+$/;
        croak "$t is not Test::Class or integer"
            unless _isa_class( __PACKAGE__, $t );
        if (my $reason = $t->SKIP_CLASS) {
            _show_header($t, @tests);
            $Builder->skip( $reason ) unless $reason eq "1";
        } else {
            $t = $t->new unless ref($t);
            my @test_methods    = _get_methods($t, TEST);
            if ( @test_methods ) {
                foreach my $method (_get_methods($t, STARTUP)) {
                    _show_header($t, @tests) unless _has_no_tests($t, $method);
                    my $method_passed = _run_method($t, $method, \@tests);
                    $all_passed = 0 unless $method_passed;
                    next TEST_OBJECT unless $method_passed;
                };
                my $class = ref($t);
                my @setup           = _get_methods($t, SETUP);
                my @teardown        = _get_methods($t, TEARDOWN);
                foreach my $test ( @test_methods ) { 
                    local $Current_method = $test;
                    $Builder->diag("\n$class->$test") if $ENV{TEST_VERBOSE};
                    my @methods_to_run = (@setup, $test, @teardown);
                    while ( my $method = shift @methods_to_run ) {
                        _show_header($t, @tests) unless _has_no_tests($t, $method);
                        $all_passed = 0 unless _run_method($t, $method, \@tests);
                        if ( _threw_exception( $t, $method ) ) {
                            my $num_to_skip = _total_num_tests($t, @methods_to_run);
                            $Builder->skip( "$method died" ) for ( 1 .. $num_to_skip );
                            last;
                        };
                    };
                };
                foreach my $method (_get_methods($t, SHUTDOWN)) {
                    _show_header($t, @tests) unless _has_no_tests($t, $method);
                    $all_passed = 0 unless _run_method($t, $method, \@tests);
                }
            }
            
        }
    }
    return($all_passed);
};

sub _find_calling_test_class {
    my $level = 0;
    while (my $class = caller(++$level)) {
        next if $class eq __PACKAGE__;
        return $class if _isa_class( __PACKAGE__, $class );
    }; 
    return(undef);
};

sub num_method_tests {
    my ($self, $method, $n) = @_;
    my $class = _find_calling_test_class( $self )
        or croak "not called in a Test::Class";
    my $info = _method_info($self, $class, $method)
        or croak "$method is not a test method of class $class";
    $info->num_tests($n) if defined($n);
    return( $info->num_tests );
};

sub num_tests {
    my $self = shift;
    croak "num_tests need to be called within a test method"
            unless defined $Current_method;
    return( $self->num_method_tests( $Current_method, @_ ) );
};

sub BAILOUT {
    my ($self, $reason) = @_;
    $Builder->BAILOUT($reason);
};

sub _last_test_if_exiting_immediately {
    $Builder->expected_tests || $Builder->current_test+1
};

sub FAIL_ALL {
    my ($self, $reason) = @_;
    my $last_test = _last_test_if_exiting_immediately();
    $Builder->expected_tests( $last_test ) unless $Builder->has_plan;
    $Builder->ok(0, $reason) until $Builder->current_test >= $last_test;
    my $num_failed = $Builder->can("history")
      ? $Builder->history->fail_count : grep( !$_, $Builder->summary );
    exit( $num_failed < 254 ? $num_failed : 254 );
};

sub SKIP_ALL {  
    my ($self, $reason) = @_;
    $Builder->skip_all( $reason ) unless $Builder->has_plan;
    my $last_test = _last_test_if_exiting_immediately();
    $Builder->skip( $reason ) 
        until $Builder->current_test >= $last_test;
    exit(0);
}

sub add_filter {
    my ( $class, $cb ) = @_;

    if ( not ref $cb eq 'CODE' ) {
        croak "Filter isn't a code-ref"
    }

    push @Filters, $cb;
}

1;

__END__

=head1 NAME

Test::Class - Easily create test classes in an xUnit/JUnit style

=head1 SYNOPSIS

  package Example::Test;
  use base qw(Test::Class);
  use Test::More;

  # setup methods are run before every test method. 
  sub make_fixture : Test(setup) {
      my $array = [1, 2];
      shift->{test_array} = $array;
  };

  # a test method that runs 1 test
  sub test_push : Test {
      my $array = shift->{test_array};
      push @$array, 3;
      is_deeply($array, [1, 2, 3], 'push worked');
  };

  # a test method that runs 4 tests
  sub test_pop : Test(4) {
      my $array = shift->{test_array};
      is(pop @$array, 2, 'pop = 2');
      is(pop @$array, 1, 'pop = 1');
      is_deeply($array, [], 'array empty');
      is(pop @$array, undef, 'pop = undef');
  };

  # teardown methods are run after every test method.
  sub teardown : Test(teardown) {
      my $array = shift->{test_array};
      diag("array = (@$array) after test(s)");
  };

later in a nearby .t file

  #! /usr/bin/perl
  use Example::Test;

  # run all the test methods in Example::Test
  Test::Class->runtests;

Outputs:

  1..5
  ok 1 - pop = 2
  ok 2 - pop = 1
  ok 3 - array empty
  ok 4 - pop = undef
  # array = () after test(s)
  ok 5 - push worked
  # array = (1 2 3) after test(s)


=head1 DESCRIPTION

Test::Class provides a simple way of creating classes and objects to test your code in an xUnit style. 

Built using L<Test::Builder>, it was designed to work with other Test::Builder based modules (L<Test::More>, L<Test::Differences>, L<Test::Exception>, etc.).

I<Note:> This module will make more sense, if you are already familiar with the "standard" mechanisms for testing perl code. Those unfamiliar with L<Test::Harness>, L<Test::Simple>, L<Test::More> and friends should go take a look at them now. L<Test::Tutorial> is a good starting point.


=head1 INTRODUCTION

=head2 A brief history lesson

In 1994 Kent Beck wrote a testing framework for Smalltalk called SUnit. It was popular. You can read a copy of his original paper at L<http://www.xprogramming.com/testfram.htm>.

Later Kent Beck and Erich Gamma created JUnit for testing Java L<http://www.junit.org/>. It was popular too.

Now there are xUnit frameworks for every language from Ada to XSLT. You can find a list at L<http://www.xprogramming.com/software.htm>.

While xUnit frameworks are traditionally associated with unit testing they are also useful in the creation of functional/acceptance tests.

Test::Class is (yet another) implementation of xUnit style testing in Perl. 


=head2 Why you should use Test::Class

Test::Class attempts to provide simple xUnit testing that integrates simply with the standard perl *.t style of testing. In particular:

=over 4

=item *

All the advantages of xUnit testing. You can easily create test fixtures and isolate tests. It provides a framework that should be familiar to people who have used other xUnit style test systems.


=item *

It is built with L<Test::Builder> and should co-exist happily with all other Test::Builder based modules. This makes using test classes in *.t scripts, and refactoring normal tests into test classes, much simpler because:

=over 4

=item *

You do not have to learn a new set of new test APIs and can continue using ok(), like(), etc. from L<Test::More> and friends. 

=item *

Skipping tests and todo tests are supported. 

=item *

You can have normal tests and Test::Class classes co-existing in the same *.t script. You don't have to re-write an entire script, but can use test classes as and when it proves useful.

=back

=item *

You can easily package your tests as classes/modules, rather than *.t scripts. This simplifies reuse, documentation and distribution, encourages refactoring, and allows tests to be extended by inheritance.

=item *

You can have multiple setup/teardown methods. For example have one teardown method to clean up resources and another to check that class invariants still hold.

=item *

It can make running tests faster. Once you have refactored your *.t scripts into classes they can be easily run from a single script. This gains you the (often considerable) start up time that each separate *.t script takes.

=back


=head2 Why you should I<not> use Test::Class

=over 4

=item *

If your *.t scripts are working fine then don't bother with Test::Class. For simple test suites it is almost certainly overkill. Don't start thinking about using Test::Class until issues like duplicate code in your test scripts start to annoy.

=item *

If you are distributing your code it is yet another module that the user has to have to run your tests (unless you distribute it with your test suite of course).

=item *

If you are used to the TestCase/Suite/Runner class structure used by JUnit and similar testing frameworks you may find Test::Unit more familiar (but try reading L</"HELP FOR CONFUSED JUNIT USERS"> before you give up).

=back


=head1 TEST CLASSES

A test class is just a class that inherits from Test::Class. Defining a test class is as simple as doing:

  package Example::Test;
  use base qw(Test::Class);

Since Test::Class does not provide its own test functions, but uses those provided by L<Test::More> and friends, you will nearly always also want to have:

  use Test::More;

to import the test functions into your test class.

=head1 METHOD TYPES

There are three different types of method you can define using Test::Class.

=head2 1) Test methods

You define test methods using the L<Test|/"Test"> attribute. For example:

  package Example::Test;
  use base qw(Test::Class);
  use Test::More;

  sub subtraction : Test {
      is( 2-1, 1, 'subtraction works );
  };

This declares the C<subtraction> method as a test method that runs one test. 

If your test method runs more than one test, you should put the number of tests in brackets like this:

  sub addition : Test(2) {
      is(10 + 20, 30, 'addition works');
      is(20 + 10, 30, '  both ways');
  };

If you don't know the number of tests at compile time you can use C<no_plan> like this.

  sub check_class : Test(no_plan) {
      my $objects = shift->{objects};
      isa_ok($_, "Object") foreach @$objects;
  };

or use the :Tests attribute, which acts just like C<:Test> but defaults to C<no_plan> if no number is given:

  sub check_class : Tests {
      my $objects = shift->{objects};
      isa_ok($_, "Object") foreach @$objects;
  };


=head2 2) Setup and teardown methods

Setup and teardown methods are run before and after every test. For example:

  sub before : Test(setup)    { diag("running before test") };
  sub after  : Test(teardown) { diag("running after test") };

You can use setup and teardown methods to create common objects used by all of your test methods (a test I<fixture>) and store them in your Test::Class object, treating it as a hash. For example:

  sub pig : Test(setup) {
      my $self = shift;
      $self->{test_pig} = Pig->new;
  };

  sub born_hungry : Test {
      my $pig = shift->{test_pig};
      is($pig->hungry, 'pigs are born hungry');
  };

  sub eats : Test(3) {
      my $pig = shift->{test_pig};
      ok(  $pig->feed,   'pig fed okay');
      ok(! $pig->hungry, 'fed pig not hungry');
      ok(! $pig->feed,   'cannot feed full pig');
  };

You can also declare setup and teardown methods as running tests. For example you could check that the test pig survives each test method by doing:

  sub pig_alive : Test(teardown => 1) {
      my $pig = shift->{test_pig};
      ok($pig->alive, 'pig survived tests' );
  };


=head2 3) Startup and shutdown methods

Startup and shutdown methods are like setup and teardown methods for the whole test class. All the startup methods are run once when you start running a test class. All the shutdown methods are run once just before a test class stops running.

You can use these to create and destroy expensive objects that you don't want to have to create and destroy for every test - a database connection for example:

  sub db_connect : Test(startup) {
      shift->{dbi} = DBI->connect;
  };

  sub db_disconnect : Test(shutdown) {
      shift->{dbi}->disconnect;
  };

Just like setup and teardown methods you can pass an optional number of tests to startup and shutdown methods. For example:

  sub example : Test(startup => 1) {
      ok(1, 'a startup method with one test');
  };
  
If you want to run an unknown number of tests within your startup method, you need to say e.g.

  sub example : Test(startup => no_plan) {
     ok(1, q{The first of many tests that don't want to have to count});
     ...
  }

as the : Tests attribute behaves exactly like : Test in this context.

If a startup method has a failing test or throws an exception then all other tests for the current test object are ignored. 

=head1 RUNNING TESTS

You run test methods with L<runtests()|"runtests">. Doing:

  Test::Class->runtests

runs all of the test methods in every loaded test class. This allows you to easily load multiple test classes in a *.t file and run them all.

  #! /usr/bin/perl
  
  # load all the test classes I want to run
  use Foo::Test;
  use Foo::Bar::Test;
  use Foo::Fribble::Test;
  use Foo::Ni::Test;
  
  # and run them all
  Test::Class->runtests;
  
You can use L<Test::Class::Load> to automatically load all the test classes in a given set of directories.

If you need finer control you can create individual test objects with L<new()|"new">. For example to just run the tests in the test class C<Foo::Bar::Test> you can do:

  Example::Test->new->runtests

You can also pass L<runtests()|/"runtests"> a list of test objects to run. For example:

  my $o1 = Example::Test->new;
  my $o2 = Another::Test->new;
  # runs all the tests in $o1 and $o2
  $o1->runtests($o2);

Since, by definition, the base Test::Class has no tests you could also have written:

  my $o1 = Example::Test->new;
  my $o2 = Another::Test->new;
  Test::Class->runtests($o1, $o2);

If you pass L<runtests()|/"runtests"> class names it will automatically create test objects for you, so the above can be written more compactly as:

  Test::Class->runtests(qw( Example::Test Another::Test ))

In all of the above examples L<runtests()|/"runtests"> will look at the number of tests both test classes run and output an appropriate test header for L<Test::Harness> automatically.

What happens if you run test classes and normal tests in the same script? For example:

  Test::Class->runtests;
  ok(Example->new->foo, 'a test not in the test class');
  ok(Example->new->bar, 'ditto');

L<Test::Harness> will complain that it saw more tests than it expected since the test header output by L<runtests()|/"runtests"> will not include the two normal tests.

To overcome this problem you can pass an integer value to L<runtests()|/"runtests">. This is added to the total number of tests in the test header. So the problematic example can be rewritten as follows:

  Test::Class->runtests(+2);
  ok(Example->new->foo, 'a test not in the test class');
  ok(Example->new->bar, 'ditto');

If you prefer to write your test plan explicitly you can use L<expected_tests()|/"expected_tests"> to find out the number of tests a class/object is expected to run.

Since L<runtests()|/"runtests"> will not output a test plan if one has already been set the previous example can be written as:

  plan tests => Test::Class->expected_tests(+2);
  Test::Class->runtests;
  ok(Example->new->foo, 'a test not in the test class');
  ok(Example->new->bar, 'ditto');

I<Remember:> Test objects are just normal perl objects. Test classes are just normal perl classes. Setup, test and teardown methods are just normal methods. You are completely free to have other methods in your class that are called from your test methods, or have object specific C<new> and C<DESTROY> methods. 

In particular you can override the new() method to pass parameters to your test object, or re-define the number of tests a method will run. See L<num_method_tests()|/"num_method_tests"> for an example. 


=head1 TEST DESCRIPTIONS

The test functions you import from L<Test::More> and other L<Test::Builder> based modules usually take an optional third argument that specifies the test description, for example:

  is $something, $something_else, 'a description of my test';
    
If you do not supply a test description, and the test function does not supply its own default, then Test::Class will use the name of the currently running test method, replacing all "_" characters with spaces so:

  sub one_plus_one_is_two : Test {
      is 1+1, 2;
  }
  
will result in:

  ok 1 - one plus one is two


=head1 RUNNING ORDER OF METHODS

Methods of each type are run in the following order:

=over 4

=item 1.

All of the startup methods in alphabetical order

=item 2.

For each test method, in alphabetical order:

=over 2

=item *

All of the setup methods in alphabetical order

=item *

The test method.

=item *

All of the teardown methods in alphabetical order

=back

=item 3.

All of the shutdown methods in alphabetical order.

=back

Most of the time you should not care what order tests are run in, but it can occasionally be useful to force some test methods to be run early. For example:

  sub _check_new {
      my $self = shift;
      isa_ok(Object->new, "Object") or $self->BAILOUT('new fails!');
  };

The leading C<_> will force the above method to run first - allowing the entire suite to be aborted before any other test methods run.


=head1 HANDLING EXCEPTIONS

If a startup, setup, test, teardown or shutdown method dies then L<runtests()|/"runtests"> will catch the exception and fail any remaining test. For example:

  sub test_object : Test(2) {
      my $object = Object->new;
      isa_ok( $object, "Object" ) or die "could not create object\n";
      ok( $object->open, "open worked" );
  };

will produce the following if the first test failed:

  not ok 1 - The object isa Object
  #   Failed test 'The object isa Object'
  #   at /Users/adrianh/Desktop/foo.pl line 14.
  #   (in MyTest->test_object)
  #     The object isn't defined
  not ok 2 - test_object died (could not create object)
  #   Failed test 'test_object died (could not create object)'
  #   at /Users/adrianh/Desktop/foo.pl line 19.
  #   (in MyTest->test_object)

This can considerably simplify testing code that throws exceptions. 

Rather than having to explicitly check that the code exited normally (e.g. with L<Test::Exception/"lives_ok">) the test will fail automatically - without aborting the other test methods. For example contrast:

  use Test::Exception;

  my $file;
  lives_ok { $file = read_file('test.txt') } 'file read';
  is($file, "content", 'test file read');

with:

  sub read_file : Test {
      is(read_file('test.txt'), "content", 'test file read');
  };
  
If more than one test remains after an exception then the first one is failed, and the remaining ones are skipped.

If the setup method of a test method dies, then all of the remaining setup and shutdown methods are also skipped.

Since startup methods will usually be creating state needed by all the other test methods an exception within a startup method will prevent all other test methods of that class running.


=head1 RETURNING EARLY

If a test method returns before it has run all of its tests, by default the missing tests are deemed to have been skipped; see L<"Skipped Tests"> for more information.

However, if the class's C<fail_if_returned_early> method returns true, then the missing tests will be deemed to have failed.  For example,

  package MyClass;
  use base 'Test::Class';
  sub fail_if_returned_early { 1 }

  sub oops : Tests(8) {
    for (my $n=1; $n*$n<50; ++$n) {
      ok 1, "$n squared is less than fifty";
    }
  }


=head1 SKIPPED TESTS

You can skip the rest of the tests in a method by returning from the method before all the test have finished running (but see L<"Returning Early"> for how to change this). The value returned is used as the reason for the tests being skipped.

This makes managing tests that can be skipped for multiple reasons very simple. For example:

  sub flying_pigs : Test(5) {
      my $pig = Pig->new;
      isa_ok($pig, 'Pig')           or return("cannot breed pigs")
      can_ok($pig, 'takeoff')       or return("pigs don't fly here");
      ok($pig->takeoff, 'takeoff')  or return("takeoff failed");
      ok( $pig->altitude > 0, 'Pig is airborne' );
      ok( $pig->airspeed > 0, '  and moving'    );
  };

If you run this test in an environment where C<Pig-E<gt>new> worked and the takeoff method existed, but failed when ran, you would get:

  ok 1 - The object isa Pig
  ok 2 - can takeoff
  not ok 3 - takeoff
  ok 4 # skip takeoff failed
  ok 5 # skip takeoff failed

You can also skip tests just as you do in Test::More or Test::Builder - see L<Test::More/"Conditional tests"> for more information. 

I<Note:> if you want to skip tests in a method with C<no_plan> tests then you have to explicitly skip the tests in the method - since Test::Class cannot determine how many tests (if any) should be skipped:

  sub test_objects : Tests {
      my $self = shift;
      my $objects = $self->{objects};
      if (@$objects) {
          isa_ok($_, "Object") foreach (@$objects);
      } else {
          $self->builder->skip("no objects to test");
      };
  };

Another way of overcoming this problem is to explicitly set the number of tests for the method at run time using L<num_method_tests()|/"num_method_tests"> or L<"num_tests">.

You can make a test class skip all of its tests by setting L<SKIP_CLASS()|SKIP_CLASS> before L<runtests()|"runtests"> is called.

=head1 TO DO TESTS

You can create todo tests just as you do in L<Test::More> and L<Test::Builder> using the C<$TODO> variable. For example:

  sub live_test : Test  {
      local $TODO = "live currently unimplemented";
      ok(Object->live, "object live");
  };

See L<Test::Harness/"Todo tests"> for more information.


=head1 EXTENDING TEST CLASSES BY INHERITANCE

You can extend test methods by inheritance in the usual way. For example consider the following test class for a C<Pig> object.

  package Pig::Test;
  use base qw(Test::Class);
  use Test::More;

  sub testing_class { "Pig" };
  sub new_args { (-age => 3) };

  sub setup : Test(setup) {
      my $self = shift;
      my $class = $self->testing_class;
      my @args = $self->new_args;
      $self->{pig} = $class->new( @args );
  };

  sub _creation : Test {
      my $self = shift;
      isa_ok($self->{pig}, $self->testing_class) 
              or $self->FAIL_ALL('Pig->new failed');
  };

  sub check_fields : Test {
      my $pig = shift->{pig};
      is($pig->age, 3, "age accessed");
  };

Next consider C<NamedPig> a subclass of C<Pig> where you can give your pig a name.

We want to make sure that all the tests for the C<Pig> object still work for C<NamedPig>. We can do this by subclassing C<Pig::Test> and overriding the C<testing_class> and C<new_args> methods.

  package NamedPig::Test;
  use base qw(Pig::Test);
  use Test::More;

  sub testing_class { "NamedPig" };
  sub new_args { (shift->SUPER::new_args, -name => 'Porky') };

Now we need to test the name method. We could write another test method, but we also have the option of extending the existing C<check_fields> method.

  sub check_fields : Test(2) {
      my $self = shift;
      $self->SUPER::check_fields;   
      is($self->{pig}->name, 'Porky', 'name accessed');
  };

While the above works, the total number of tests for the method is dependent on the number of tests in its C<SUPER::check_fields>. If we add a test to C<Pig::Test-E<gt>check_fields> we will also have to update the number of tests of C<NamedPig::test-E<gt>check_fields>.

Test::Class allows us to state explicitly that we are adding tests to an existing method by using the C<+> prefix. Since we are adding a single test to C<check_fields> it can be rewritten as:

  sub check_fields : Test(+1) {
      my $self = shift;
      $self->SUPER::check_fields;
      is($self->{pig}->name, 'Porky', 'name accessed');
  };

With the above definition you can add tests to C<check_fields> in C<Pig::Test> without affecting C<NamedPig::Test>.


=head1 RUNNING INDIVIDUAL TESTS

B<NOTE:> The exact mechanism for running individual tests is likely to change in the future. 

Sometimes you just want to run a single test.  Commenting out other tests or writing code to skip them can be a hassle, so you can specify the C<TEST_METHOD> environment variable.  The value is expected to be a valid regular expression and, if present, only runs test methods whose names match the regular expression.  Startup, setup, teardown and shutdown tests will still be run.

One easy way of doing this is by specifying the environment variable I<before> the C<runtests> method is called.

Running a test named C<customer_profile>:

 #! /usr/bin/perl
 use Example::Test;
      
 $ENV{TEST_METHOD} = 'customer_profile';
 Test::Class->runtests;

Running all tests with C<customer> in their name:

 #! /usr/bin/perl
 use Example::Test;
      
 $ENV{TEST_METHOD} = '.*customer.*';
 Test::Class->runtests;

If you specify an invalid regular expression, your tests will not be run:

 #! /usr/bin/perl
 use Example::Test;
      
 $ENV{TEST_METHOD} = 'C++';
 Test::Class->runtests;

And when you run it:

 TEST_METHOD (C++) is not a valid regular expression: Search pattern \
 not terminated at (eval 17) line 1.


=head1 ORGANISING YOUR TEST CLASSES

You can, of course, organise your test modules as you wish. My personal preferences is:

=over 4

=item *

Name test classes with a suffix of C<::Test> so the test class for the C<Foo::Bar> module would be C<Foo::Bar::Test>.

=item *

Place all test classes in F<t/lib>.

=back

The L<Test::Class::Load> provides a simple mechanism for easily loading all of the test classes in a given set of directories.


=head1 A NOTE ON LOADING TEST CLASSES

Due to its use of subroutine attributes Test::Class based modules must be loaded at compile rather than run time. This is because the :Test attribute is applied by a CHECK block.

This can be problematic if you want to dynamically load Test::Class modules. Basically while:

  require $some_test_class;
  
will break, doing:

  BEGIN { require $some_test_class };
  
will work just fine. For more information on CHECK blocks see L<perlmod/"BEGIN, CHECK, INIT and END">.

If you still can't arrange for your classes to be loaded at runtime, you could use an alternative mechanism for adding your tests:

  # sub test_something : Test(3) {...}
  # becomes
  sub test_something {...}
  __PACKAGE__->add_testinfo('test_something', test => 3);

See the L<add_testinfo|/"add_testinfo"> method for more details.

=head1 GENERAL FILTERING OF TESTS

The use of $ENV{TEST_METHOD} to run just a subset of tests is useful, but
sometimes it doesn't give the level of granularity that you desire.  Another
feature of this class is the ability to do filtering on other static criteria.
In order to permit this, a generic filtering method is supported.  This can
be used by specifying coderefs to the 'add_filter' method of this class.

In determining which tests should be run, all filters that have previously
been specified via the add_filter method will be run in-turn for each normal
test method.  If B<any> of these filters return a false value, the method will
not be executed, or included in the number of tests.  Note that filters will
only be run for normal test methods, they are ignored for startup, shutdown,
setup, and teardown test methods.

Note that test filters are global, and will affect all tests in all classes,
not just the one that they were defined in.

An example of this mechanism that mostly simulates the use of TEST_METHOD
above is:

 package MyTests;

 use Test::More;

 use base qw( Test::Class );

 my $MYTEST_METHOD = qr/^t_not_filtered$/;

 my $filter = sub {
    my ( $test_class, $test_method ) = @_;

    return $test_method =~ $MYTEST_METHOD;
 };
 Test::Class->add_filter( $filter );

 sub t_filtered : Test( 1 ) {
    fail( "filtered test run" );
 }

 sub t_not_filtered : Test( 1 ) {
    pass( "unfiltered test run" );
 }

=head1 METHODS

=head2 Creating and running tests

=over 4

=item B<Test>

  # test methods
  sub method_name : Test { ... };
  sub method_name : Test(N) { ... };

  # setup methods
  sub method_name : Test(setup) { ... };
  sub method_name : Test(setup => N) { ... };

  # teardown methods
  sub method_name : Test(teardown) { ... };
  sub method_name : Test(teardown => N) { ... };

  # startup methods
  sub method_name : Test(startup) { ... };
  sub method_name : Test(startup => N) { ... };

  # shutdown methods
  sub method_name : Test(shutdown) { ... };
  sub method_name : Test(shutdown => N) { ... };

Marks a startup, setup, test, teardown or shutdown method. See L<runtests()|/"runtests"> for information on how to run methods declared with the C<Test> attribute.

N specifies the number of tests the method runs. 

=over 4

=item *

If N is an integer then the method should run exactly N tests.

=item *

If N is an integer with a C<+> prefix then the method is expected to call its C<SUPER::> method and extend it by running N additional tests.

=item *

If N is the string C<no_plan> then the method can run an arbitrary number of tests.

=back

If N is not specified it defaults to C<1> for test methods, and C<0> for startup, setup, teardown and shutdown methods. 

You can change the number of tests that a method runs using L<num_method_tests()|/"num_method_tests"> or L<num_tests()|/"num_tests">.


=item B<Tests>

  sub method_name : Tests { ... };
  sub method_name : Tests(N) { ... };

Acts just like the C<:Test> attribute, except that if the number of tests is not specified it defaults to C<no_plan>. So the following are equivalent:

  sub silly1 :Test( no_plan ) { ok(1) foreach (1 .. rand 5) };
  sub silly2 :Tests           { ok(1) foreach (1 .. rand 5) };


=item B<new>

  $Tests = CLASS->new(KEY => VAL ...)
  $Tests2 = $Tests->new(KEY => VAL ...)

Creates a new test object (blessed hashref) containing the specified key/value pairs. 

If called as an object method the existing object's key/value pairs are copied into the new object. Any key/value pairs passed to C<new> override those in the original object if duplicates occur.

Since the test object is passed to every test method as it runs it is a convenient place to store test fixtures. For example:

  sub make_fixture : Test(setup) {
      my $self = shift;
      $self->{object} = Object->new();
      $self->{dbh} = Mock::DBI->new(-type => normal);
  };

  sub test_open : Test {
      my $self = shift;
      my ($o, $dbh) = ($self->{object}, $self->{dbh});
      ok($o->open($dbh), "opened ok");
  };

See L<num_method_tests()|/"num_method_tests"> for an example of overriding C<new>.


=item B<expected_tests>

  $n = $Tests->expected_tests
  $n = CLASS->expected_tests
  $n = $Tests->expected_tests(TEST, ...)
  $n = CLASS->expected_tests(TEST, ...)

Returns the total number of tests that L<runtests()|/"runtests"> will run on the specified class/object. This includes tests run by any setup and teardown methods.

Will return C<no_plan> if the exact number of tests is undetermined (i.e. if any setup, test or teardown method has an undetermined number of tests).

The C<expected_tests> of an object after L<runtests()|/"runtests"> has been executed will include any run time changes to the expected number of tests made by L<num_tests()|/"num_tests"> or L<num_method_tests()|/"num_method_tests">.

C<expected_tests> can also take an optional list of test objects, test classes and integers. In this case the result is the total number of expected tests for all the test/object classes (including the one the method was applied to) plus any integer values.

C<expected_tests> is useful when you're integrating one or more test classes into a more traditional test script, for example:

  use Test::More;
  use My::Test::Class;

  plan tests => My::Test::Class->expected_tests(+2);

  ok(whatever, 'a test');
  ok(whatever, 'another test');
  My::Test::Class->runtests;



=item B<runtests>

  $allok = $Tests->runtests
  $allok = CLASS->runtests
  $allok = $Tests->runtests(TEST, ...)
  $allok = CLASS->runtests(TEST, ...)

C<runtests> is used to run test classes. At its most basic doing:

  $test->runtests
  
will run the test methods of the test object $test, unless C<< $test->SKIP_CLASS >> returns a true value. 

Unless you have already specified a test plan using Test::Builder (or Test::More, et al) C<runtests> will set the test plan just before the first method that runs a test is executed. 

If the environment variable C<TEST_VERBOSE> is set C<runtests> will display the name of each test method before it runs like this:

  # My::Test::Class->my_test
  ok 1 - fribble
  # My::Test::Class->another_test
  ok 2 - bar

Just like L<expected_tests()|/"expected_tests">, C<runtests> can take an optional list of test object/classes and integers. All of the test object/classes are run. Any integers are added to the total number of tests shown in the test header output by C<runtests>. 

For example, you can run all the tests in test classes A, B and C, plus one additional normal test by doing:

    Test::Class->runtests(qw(A B C), +1);
    ok(1==1, 'non class test');

Finally, if you call C<runtests> on a test class without any arguments it will run all of the test methods of that class, and all subclasses of that class. For example:

  #! /usr/bin/perl
  # Test all the Foo stuff
  
  use Foo::Test;
  use Foo::Bar::Test;
  use Foo::Ni::Test;
  
  # run all the Foo*Test modules we just loaded
  Test::Class->runtests;
    

=item B<SKIP_CLASS>

  $reason = CLASS->SKIP_CLASS;
  CLASS->SKIP_CLASS( $reason );

Determines whether the test class CLASS should run it's tests. If SKIP_CLASS returns a true value then  L<runtests()|/"runtests"> will not run any of the test methods in CLASS.

You can override the default on a class-by-class basis by supplying a new value to SKIP_CLASS. For example if you have an abstract base class that should not run just add the following to your module:

  My::Abstract::Test->SKIP_CLASS( 1 );
  
This will not affect any sub-classes of C<My::Abstract::Test> which will run as normal.

If the true value returned by SKIP_CLASS is anything other than "1" then a skip test is output using this value as the skip message. For example:

  My::Postgres::Test->SKIP_CLASS(
      $ENV{POSTGRES_HOME} ? 0 : '$POSTGRES_HOME needs to be set'
  );

will output something like this if C<POSTGRES_HOME> is not set

    ... other tests ...
    ok 123 # skip My::Postgres::Test  - $POSTGRES_HOME needs to be set
    ... more tests ...

You can also override SKIP_CLASS for a class hierarchy. For example, to prevent any subclasses of My::Postgres::Test running we could override SKIP_CLASS like this:

  sub My::Postgres::Test::SKIP_CLASS {
      $ENV{POSTGRES_HOME} ? 0 : '$POSTGRES_HOME needs to be set'
  };
  
=back

=head2 Fetching and setting a method's test number


=over 4

=item B<num_method_tests>

  $n = $Tests->num_method_tests($method_name)
  $Tests->num_method_tests($method_name, $n)
  $n = CLASS->num_method_tests($method_name)
  CLASS->num_method_tests($method_name, $n)

Fetch or set the number of tests that the named method is expected to run.

If the method has an undetermined number of tests then $n should be the string C<no_plan>.

If the method is extending the number of tests run by the method in a superclass then $n should have a C<+> prefix.

When called as a class method any change to the expected number of tests applies to all future test objects. Existing test objects are unaffected. 

When called as an object method any change to the expected number of tests applies to that object alone.

C<num_method_tests> is useful when you need to set the expected number of tests at object creation time, rather than at compile time.

For example, the following test class will run a different number of tests depending on the number of objects supplied.

  package Object::Test; 
  use base qw(Test::Class);
  use Test::More;

  sub new {
      my $class = shift;
      my $self = $class->SUPER::new(@_);
      my $num_objects = @{$self->{objects}};
      $self->num_method_tests('test_objects', $num_objects);
      return($self);
  };

  sub test_objects : Tests {
    my $self = shift;
    ok($_->open, "opened $_") foreach @{$self->{objects}};
  };
  ...
  # This runs two tests
  Object::Test->new(objects => [$o1, $o2]);

The advantage of setting the number of tests at object creation time, rather than using a test method without a plan, is that the number of expected tests can be determined before testing begins. This allows better diagnostics from L<runtests()|/"runtests">, L<Test::Builder> and L<Test::Harness>.

C<num_method_tests> is a protected method and can only be called by subclasses of Test::Class. It fetches or sets the expected number of tests for the methods of the class it was I<called in>, not the methods of the object/class it was I<applied to>. This allows test classes that use C<num_method_tests> to be subclassed easily.

For example, consider the creation of a subclass of Object::Test that ensures that all the opened objects are read-only:

  package Special::Object::Test;
  use base qw(Object::Test);
  use Test::More;

  sub test_objects : Test(+1) {
      my $self = shift;
      $self->SUPER::test_objects;
      my @bad_objects = grep {! $_->read_only} (@{$self->{objects}});
      ok(@bad_objects == 0, "all objects read only");
  };
  ...
  # This runs three tests
  Special::Object::Test->new(objects => [$o1, $o2]);

Since the call to C<num_method_tests> in Object::Test only affects the C<test_objects> of Object::Test, the above works as you would expect.


=item B<num_tests>

  $n = $Tests->num_tests
  $Tests->num_tests($n)
  $n = CLASS->num_tests
  CLASS->num_tests($n)

Set or return the number of expected tests associated with the currently running test method. This is the same as calling L<num_method_tests()|/"num_method_tests"> with a method name of L<current_method()|/"current_method">.

For example:

  sub txt_files_readable : Tests {
      my $self = shift;
      my @files = <*.txt>;
      $self->num_tests(scalar(@files));
      ok(-r $_, "$_ readable") foreach (@files);
  };

Setting the number of expected tests at run time, rather than just having a C<no_plan> test method, allows L<runtests()|/"runtests"> to display appropriate diagnostic messages if the method runs a different number of tests.

=back


=head2 Support methods

=over 4

=item B<builder>

  $Tests->builder

Returns the underlying L<Test::Builder> object that Test::Class uses. For example:

  sub test_close : Test {
      my $self = shift;
      my ($o, $dbh) = ($self->{object}, $self->{dbh});
      $self->builder->ok($o->close($dbh), "closed ok");
  };

=item B<current_method>

  $method_name = $Tests->current_method
  $method_name = CLASS->current_method

Returns the name of the test method currently being executed by L<runtests()|/"runtests">, or C<undef> if L<runtests()|/"runtests"> has not been called. 

The method name is also available in the setup and teardown methods that run before and after the test method. This can be useful in producing diagnostic messages, for example:

  sub test_invarient : Test(teardown => 1) {
      my $self = shift;
      my $m = $self->current_method;
      ok($self->invarient_ok, "class okay after $m");
  };



=item B<BAILOUT>

  $Tests->BAILOUT($reason)
  CLASS->BAILOUT($reason)

Things are going so badly all testing should terminate, including running any additional test scripts invoked by L<Test::Harness>. This is exactly the same as doing:

  $self->builder->BAILOUT

See L<Test::Builder/"BAILOUT"> for details. Any teardown and shutdown methods are I<not> run.


=item B<FAIL_ALL>

  $Tests->FAIL_ALL($reason)
  CLASS->FAIL_ALL($reason)

Things are going so badly all the remaining tests in the current script should fail. Exits immediately with the number of tests failed, or C<254> if more than 254 tests were run. Any teardown methods are I<not> run.

This does not affect the running of any other test scripts invoked by L<Test::Harness>.

For example, if all your tests rely on the ability to create objects then you might want something like this as an early test:

  sub _test_new : Test(3) {
      my $self = shift;
      isa_ok(Object->new, "Object") 
          || $self->FAIL_ALL('cannot create Objects');
      ...
  };



=item B<SKIP_ALL>

  $Tests->SKIP_ALL($reason)
  CLASS->SKIP_ALL($reason)

Things are going so badly all the remaining tests in the current script should be skipped. Exits immediately with C<0> - teardown methods are I<not> run.

This does not affect the running of any other test scripts invoked by L<Test::Harness>.

For example, if you had a test script that only applied to the darwin OS you could write:

  sub _darwin_only : Test(setup) {
      my $self = shift;
      $self->SKIP_ALL("darwin only") unless $^O eq "darwin";    
  };


=item B<add_testinfo>

  CLASS->add_testinfo($name, $type, $num_tests)

Chiefly for use by libraries like L<Test::Class::Sugar>, which can't use the C<:Test(...)> interfaces make test methods. C<add_testinfo> informs the class about a test method that has been defined without a C<Test>, C<Tests> or other attribute.

C<$name> is the name of the method, C<$type> must be one of C<startup>, C<setup>, C<test>, C<teardown> or C<shutdown>, and C<$num_tests> has the same meaning as C<N> in the description of the L<Test|/"Test"> attribute.


=item B<add_filter>

    CLASS->add_filter($filter_coderef);

Adds a filtering coderef. Each filter is passed a test class and method name and returns a boolean. All filters are applied globally in the order they were added. If any filter returns false the test method is not run or included in the number of tests.

Note that filters will only be run for normal test methods, they are ignored for startup, shutdown, setup, and teardown test methods.

See the section on the L</"GENERAL FILTERING OF TESTS"> for more information.

=item B<fail_if_returned_early>

Controls what happens if a method returns before it has run all of its tests.  It is called with no arguments in boolean context; if it returns true, then the missing tests fail, otherwise, they skip.  See L<"Returning Early"> and L<"Skipped Tests">.


=back



=head1 HELP FOR CONFUSED JUNIT USERS

This section is for people who have used JUnit (or similar) and are confused because they don't see the TestCase/Suite/Runner class framework they were expecting. Here we take each of the major classes in JUnit and compare them with their equivalent Perl testing modules.

=over 4

=item B<Class Assert>

The test assertions provided by Assert correspond to the test functions provided by the L<Test::Builder> based modules (L<Test::More>, L<Test::Exception>, L<Test::Differences>, etc.)

Unlike JUnit the test functions supplied by Test::More et al do I<not> throw exceptions on failure. They just report the failure to STDOUT where it is collected by L<Test::Harness>. This means that where you have

  sub foo : Test(2) {
      ok($foo->method1);
      ok($foo->method2);
  };

The second test I<will> run if the first one fails. You can emulate the JUnit way of doing it by throwing an explicit exception on test failure:

  sub foo : Test(2) {
      ok($foo->method1) or die "method1 failed";
      ok($foo->method2);
  };

The exception will be caught by Test::Class and the other test automatically failed.

=item B<Class TestCase>

Test::Class corresponds to TestCase in JUnit.

In Test::Class setup, test and teardown methods are marked explicitly using the L<Test|/"Test"> attribute. Since we need to know the total number of tests to provide a test plan for L<Test::Harness> we also state how many tests each method runs.

Unlike JUnit you can have multiple setup/teardown methods in a class.

=item B<Class TestSuite>

Test::Class also does the work that would be done by TestSuite in JUnit.

Since the methods are marked with attributes Test::Class knows what is and isn't a test method. This allows it to run all the test methods without having the developer create a suite manually, or use reflection to dynamically determine the test methods by name. See the L<runtests()|/"runtests"> method for more details.

The running order of the test methods is fixed in Test::Class. Methods are executed in alphabetical order.

Unlike JUnit, Test::Class currently does not allow you to run individual test methods.

=item B<Class TestRunner>

L<Test::Harness> does the work of the TestRunner in JUnit. It collects the test results (sent to STDOUT) and collates the results.

Unlike JUnit there is no distinction made by Test::Harness between errors and failures. However, it does support skipped and todo test - which JUnit does not.

If you want to write your own test runners you should look at L<Test::Harness::Straps>.

=back


=head1 OTHER MODULES FOR XUNIT TESTING IN PERL

In addition to Test::Class there are two other distributions for xUnit testing in perl. Both have a longer history than Test::Class and might be more suitable for your needs. 

I am biased since I wrote Test::Class - so please read the following with appropriate levels of scepticism. If you think I have misrepresented the modules please let me know.

=over 4

=item B<Test::SimpleUnit>

A very simple unit testing framework. If you are looking for a lightweight single module solution this might be for you.

The advantage of L<Test::SimpleUnit> is that it is simple! Just one module with a smallish API to learn. 

Of course this is also the disadvantage. 

It's not class based so you cannot create testing classes to reuse and extend.

It doesn't use L<Test::Builder> so it's difficult to extend or integrate with other testing modules. If you are already familiar with L<Test::Builder>, L<Test::More> and friends you will have to learn a new test assertion API. It does not support L<todo tests|Test::Harness/"Todo tests">.

=item B<Test::Unit>

L<Test::Unit> is a port of JUnit L<http://www.junit.org/> into perl. If you have used JUnit then the Test::Unit framework should be very familiar.

It is class based so you can easily reuse your test classes and extend by subclassing. You get a nice flexible framework you can tweak to your heart's content. If you can run Tk you also get a graphical test runner.

However, Test::Unit is not based on L<Test::Builder>. You cannot easily move Test::Builder based test functions into Test::Unit based classes. You have to learn another test assertion API. 

Test::Unit implements it's own testing framework separate from L<Test::Harness>. You can retrofit *.t scripts as unit tests, and output test results in the format that L<Test::Harness> expects, but things like L<todo tests|Test::Harness/"Todo tests"> and L<skipping tests|Test::Harness/"Skipping tests"> are not supported. 

=back


=head1 BUGS

None known at the time of writing.

If you find any bugs please let me know by e-mail at <adrianh@quietstars.com>, or report the problem with L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Test-Class>.


=head1 COMMUNITY

=head2 perl-qa

If you are interested in testing using Perl I recommend you visit L<http://qa.perl.org/> and join the excellent perl-qa mailing list. See L<http://lists.perl.org/showlist.cgi?name=perl-qa> for details on how to subscribe.

=head2 perlmonks

You can find users of Test::Class, including the module author, on  L<http://www.perlmonks.org/>. Feel free to ask questions on Test::Class there.

=head2 CPAN::Forum

The CPAN Forum is a web forum for discussing Perl's CPAN modules.   The Test::Class forum can be found at L<http://www.cpanforum.com/dist/Test-Class>.


=head1 TO DO

If you think this module should do something that it doesn't (or does something that it shouldn't) please let me know.

You can see my current to do list at L<http://adrianh.tadalist.com/lists/public/4798>, with an RSS feed of changes at L<http://adrianh.tadalist.com/lists/feed_public/4798>.


=head1 ACKNOWLEDGMENTS

This is yet another implementation of the ideas from Kent Beck's Testing Framework paper L<http://www.xprogramming.com/testfram.htm>.

Thanks to 
Adam Kennedy,
agianni,
Alexander D'Archangel,
Andrew Grangaard,
Apocalypse,
Ask Bjorn Hansen,
Chris Dolan,
Chris Williams,
Corion, 
Cosimo Streppone,
Daniel Berger,
Dave Evans,
Dave O'Neill,
David Cantrell,
David Wheeler,
Diab Jerius,
Emil Jansson, 
Gunnar Wolf,
Hai Pham,
Hynek,
imacat,
Jeff Deifik,
Jim Brandt,
Jochen Stenzel,
Johan Lindstrom, 
John West,
Jonathan R. Warden,
Joshua ben Jore,
Jost Krieger,
Ken Fox,
Kenichi Ishigaki
Lee Goddard,
Mark Morgan,
Mark Reynolds,
Mark Stosberg,
Martin Ferrari,
Mathieu Sauve-Frankel,
Matt Trout,
Matt Williamson,
Michael G Schwern, 
Murat Uenalan, 
Naveed Massjouni,
Nicholas Clark,
Ovid, 
Piers Cawley,
Rob Kinyon,
Sam Raymer,
Scott Lanning,
Sebastien Aperghis-Tramoni,
Steve Kirkup,
Stray Toaster,
Ted Carnahan,
Terrence Brannon, 
Todd W,
Tom Metro,
Tony Bowden, 
Tony Edwardson,
William McKee, 
various anonymous folk and all the fine people on perl-qa for their feedback, patches, suggestions and nagging.

This module wouldn't be possible without the excellent L<Test::Builder>. Thanks to chromatic and Michael G Schwern for creating such a useful module.


=head1 AUTHORS

Adrian Howard <adrianh@quietstars.com>, Curtis "Ovid" Poe, <ovid at cpan.org>, Mark Morgan <makk384@gmail.com>.

If you use this module, and can spare the time please let us know or rate it at L<http://cpanratings.perl.org/rate/?distribution=Test-Class>.

=head1 SEE ALSO

=over 4

=item L<Test::Class::Load>

Simple way to load "Test::Class" classes automatically.

=item L<http://del.icio.us/tag/Test::Class>

Delicious links on Test::Class.

=item Perl Testing: A Developer's Notebook by Ian Langworth and chromatic

Chapter 8 covers using Test::Class.

=item Advanced Perl Programming, second edition by Simon Cozens

Chapter 8 has a few pages on using Test::Class.

=item The Perl Journal, April 2003

Includes the article "Test-Driven Development in Perl" by Piers Cawley that uses Test::Class.

=item L<Test::Builder>

Support module for building test libraries.

=item L<Test::Simple> & L<Test::More>

Basic utilities for writing tests.

=item L<http://qa.perl.org/test-modules.html>

Overview of some of the many testing modules available on CPAN.

=item L<http://del.icio.us/tag/perl+testing>

Delicious links on perl testing.

=item L<Test::Object>

Another approach to object oriented testing.

=item L<Test::Group> and L<Test::Block>

Alternatives to grouping sets of tests together.

=back

The following modules use Test::Class as part of their test suite. You might want to look at them for usage examples:

=over 4

L<Aspect>, Bricolage (L<http://www.bricolage.cc/>), L<Class::StorageFactory>, L<CGI::Application::Search>, L<DBIx::Romani>, L<Xmldoom>, L<Object::Relational>, L<File::Random>, L<Geography::JapanesePrefectures>, L<Google::Adwords>, L<Merge::HashRef>, L<PerlBuildSystem>, L<Pixie>, L<Yahoo::Marketing>, and L<XUL-Node>

=back

The following modules are not based on L<Test::Builder>, but may be of interest as alternatives to Test::Class.

=over 4

=item L<Test::Unit>

Perl unit testing framework closely modeled on JUnit. 

=item L<Test::SimpleUnit>

A very simple unit testing framework. 

=back

=head1 LICENCE

Copyright 2002-2010 Adrian Howard, All Rights Reserved.

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.
