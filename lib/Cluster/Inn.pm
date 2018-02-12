=pod

THIS IS AN IMPLEMENTATION OF THE INN CLUSTERING ALGORITHM WITH A PRIORI RADIUS
ASSIGNMENT.

=cut

package Cluster::Inn;
use Moose;
#use lib '/home/nfulton/cs/geometry/lib';


with 'Cluster';

sub do_iteration {
    my($self , $radius) = @_;

    my $converges = 1;

    #allow each centroid to choose the best point.
    for(my $i=0;$i<scalar @{$self->centroids};$i++) {
        my $closest = $self->_find_closest_point($i);
        if(defined $closest and $closest != -1) {
            $self->_set_association($closest, $i);
            $converges = 0;
        }
    }

    my @old_centroids = @{$self->centroids};
    $self->_move_centroids($radius);

    return $converges;
}

## START FIND CLOSEST POINT ##
sub _find_closest_point {
    my($self, $centroid_index) = @_;
    my $centroid = $self->centroids->[$centroid_index];

    my $min_distance = -1;
    my $best_point   = 0;

    for(my $point_index = 0 ; $point_index < scalar @{$self->points}; $point_index++) {
        my $point    = $self->points->[$point_index];
        my $distance = $self->get_distance($point, $centroid);

        #ensure that we don't already own this point.
        my $current_centroid_index = $self->get_assoc_for_pt($point_index);
        next if $current_centroid_index == $centroid_index;

        #ensure we are only stealing when we really are the better person.
        if(defined $current_centroid_index and $current_centroid_index != -1) {
            my $current_distance = $self->get_distance($point,$self->centroids->[$current_centroid_index]);
            #print $distance . " " . $current_distance . "\n" if $distance <= $current_distance;
            next if $distance >= $current_distance;
        }

        #ensure that we're not wiping out a better choice from earlier.
        if($min_distance == -1 || $distance < $min_distance) {
            $min_distance = $distance;
            $best_point   = $point_index;
        }
    }

    if($min_distance == -1) {
        return;
    }
    return $best_point;
}


##START SET ASSOC ##
#set = update xor add.
sub _set_association {
    my $self  = shift;
    my $pt = shift;
    my $cluster = shift;
    
    my $done = 0;
    map {
        if($_->[0] == $pt) {
            $_->[1] = $cluster;
            $done = 1;
        }
    } @{$self->associations};

    if($done == 0) {
        my @new_data = ($pt, $cluster);
        push @{$self->associations} , \@new_data
    }
}

## START MOVE CENTROIDS ##
sub _move_centroids {
    my($self, $radius) = @_;
    for(my $i=0;$i<scalar @{$self->centroids};$i++) {
        $self->_move_centroid($i, $radius);
    }
}
sub _move_centroid {
    my($self, $centroid_index, $point_radius) = @_;

    my @points;
    foreach my $assoc (@{$self->associations}) {
        push @points, $assoc->[0] if($assoc->[1] == $centroid_index);
    }

    my $x = 0;
    my $y = 0;
    my $m = 0;
    foreach my $point_index (@points) {
        my $mass = $self->_get_point_mass($point_index, $point_radius);
        $x += $mass * $self->points->[$point_index]->[0];
        $y += $mass * $self->points->[$point_index]->[1];
        $m += $mass;
    }

    my @center;
    if($m != 0) {
        @center = ($x/$m , $y/$m);
    } 
    else { 
        @center = (0,0); 
    }

    $self->centroids->[$centroid_index] = \@center;
}

sub _get_point_mass {
    my($self,$point_index,$radius) = @_;
    my $point = $self->points->[$point_index];

    #total weight will = the combined differences between the radius and the
    #length to a point in that circle.
    my $total_weight = 1; 
    map {
        my $curr_pt = $_;
        die unless $curr_pt and $self->points->[$point_index];

        my $weight = $radius - $self->get_distance($self->points->[$point_index], $curr_pt);
        if($weight > 0) {
            $total_weight += $weight;
        }
    } @{$self->points};
    
    return $total_weight;
}


1;

