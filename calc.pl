use Data::Dumper;
use Math::Integral::Romberg 'integral';

#finds the centroid of a region generically.
sub _find_lamina_centroid {
    my($self, $p_ref, $top_ref, $bottom_ref, $a, $b) = @_;

    my $m_ref = sub {
        my $x = shift;
        my $ret_var = ($p_ref->($x) * ($top_ref->($x) - $bottom_ref->($x)));
        return $ret_var;
    };
    my $m = integral($m_ref , $a , $b);

    #My
    my $x_ref = sub {
        my $x = shift;
        my $ret_var = ($p_ref->($x) * $x * ($top_ref->($x) - $bottom_ref->($x)));
        print $x . " " ;
        return $ret_var;
    };
    my $x = integral($x_ref , $a , $b); 
    
    #Mx
    my $y_ref = sub {
        my $x = shift;
        my $topv = $top_ref->($x);
        $topv = $topv*$topv;
        my $bottomv = $bottom_ref->($x);
        $bottomv = $bottomv * $bottomv;
        return (($p_ref->($x)/2) * ($topv - $bottomv));
    };
    my $y = integral($y_ref , $a, $b);

    my @center = ($x/$m , $y/$m);
    return \@center;
}

sub move_center {
    my($self, $centroid_index, $centroid_radius, $point_radius) = @_;

    my $p_ref = sub {
        my $x = shift;
        my @point = ()
        return $self->get_mass_at_point(

    $self->_find_lamina_centroid(




