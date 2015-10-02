#==========================================================================
# Simple attribute class to facilitate easier test parsing..
#==========================================================================
package CESMTest;

sub new
{
    my ($class, %params) = @_;

    my $self = {
        compset => $params{'compset'} || undef,
        grid    => $params{'grid'} || undef,
        testname => $params{'testname'} || undef,
        machine  => $params{'machine'} || undef,
        compiler  => $params{'compiler'} || undef,
        testmods  => $params{'testmods'} || undef,
        comment  => $params{'comment'} || undef,
    };
    bless $self, $class;
    return $self;
}
1;
