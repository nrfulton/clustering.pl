
use strict;
use warnings;

package Cluster; 
use Moose;
use Data::Dumper;
use  Carp qw(cluck);
use Statistics::R;
use Math::Round;

#point = a place on the graph
#node  = a member of a cluster
#a node is a point bot not necessarily vice versa.

#A list of the coordinates of all nodes.
has 'points' => (
    is  => 'rw',
    isa => 'ArrayRef',
    default => sub {
        my @array;
        return \@array;
    },
);

#A list of the coordinates of centroids.
has 'centroids' => (
    is  => 'rw',
    isa => 'ArrayRef',
    default => sub {
        my @array;
        return \@array;
    },
);

#resolves a node to a centroid.
#Array of arrays w/ [0]=node index and [1]=centroid index
has 'associations' => (
    is  => 'rw',
    isa => 'ArrayRef',
    default => sub {
        my @array;
        return \@array;
    },
);

#defaults to pythagorean spaces
has 'distance_equation' => (
    #isa     => 'SubRef',
    is      => 'rw',
    default => sub {
        return sub {
            my($c1,$c2)=@_;
            cluck('err') unless ref $c1 eq "ARRAY" and ref $c2 eq "ARRAY";
            return sqrt(abs(($c2->[1] - $c1->[1]) + ($c2->[0] - $c1->[0])));
        };
    },
);
    sub get_distance {
        my $self = shift;
        return $self->distance_equation->(@_);
    }


#begin vanilla manipulators

#nodes can be added, but not changed.
sub add_point {
    my $self = shift;
    push @{$self->points} , \@_;
}

sub add_centroid {
    my $self = shift;
    push @{$self->centroids} , \@_;
}

sub add_association {
    my $self = shift;
    push @{$self->associations} , \@_;
}

sub association_at {
    my($self, $point) = @_;
    map {
        return $_->[1] if $_->[0] == $point;
    } @{$self->associations};
    return -1;
}

sub change_centroid {
    my $self  = shift;
    my $index = shift;
    $self->centroids->[$index] = \@_;
}

sub set_association {
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
        $self->add_association($pt, $cluster);
    }
}

sub get_assoc_for_pt {
    my $self = shift;
    my $pt = shift;

    map {
        return $_->[1] if $_->[0] == $pt;
    } @{$self->associations};

    return -1;
}

#begin helper algorithms

sub get_cluster_density {
    my($self, $centroid_index) = @_;

    my $count = 1; #prevent division by 0
    map {
        $count += 1 if $_->[1] == $centroid_index;
    } @{$self->associations};

    return $count;
}

sub get_mass_at_point {
    my($self, $point, $radius) = @_;
    
    my $total_mass = 1;
    map {
        my $curr_pt = $_;
        die unless $curr_pt and $point;

        my $mass = $radius - $self->get_distance($_, $point);
        if($mass > 0) {
            $total_mass += $mass;
        }
    } @{$self->points};

    return $total_mass;
}



#can pass a point index or x/y coords.
sub get_point_density {
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




#TODO this is stuff that would be moved into centroid and node classes if I
#weren't so lazy.
sub find_closest_centroid {
    my($self , $point) = @_;
    
    my $min_distance  = 0;
    my $best_centroid = 0;
    for(my $i = 0 ; $i < scalar @{$self->centroids} ; $i++) {
        my $d = $self->get_distance($point , $self->centroids->[$i]);
        if($min_distance == 0 || $d < $min_distance) {
            $min_distance  = $d;
            $best_centroid = $i;
        }
    }

    return $best_centroid;
}

sub find_closest_point {
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

=pod
    my $min_distance  = -1;
    my $best_point    = 0;
    for(my $i = 0 ; $i < scalar @{$self->points} ; $i++) {
        my $d = $self->get_distance($self->centroids->[$centroid] , $self->points->[$i]);
        if(($min_distance == -1 or $d < $min_distance) and
            $self->association_at($i) != $centroid) {
            
            #if there's a current centroid, make sure we deal w/ it properly.
            my $current_centroid = $self->get_assoc_for_pt($i);
            if( $current_centroid != -1 ) {
                #ensure all the points are assigned before we start stealing.
                next if scalar @{$self->associations} != scalar @{$self->points};
                #only steal if we really are better than the other guy.
                next if $self->get_distance($self->centroids->[$current_centroid],$self->points->[$i]) < $self->get_distance($self->centroids->[$centroid],$self->points->[$i]);
                print "centroid $centroid is stealing " . $self->points->[$i]->[0] . "," . $self->points->[$i]->[0] .
                " from $current_centroid B/C " .
                $self->get_distance($self->centroids->[$current_centroid],$self->points->[$i])
                . " > " .
                $self->get_distance($self->centroids->[$centroid],$self->points->[$i])
                . "\n";
            }

            $min_distance  = $d;
            $best_point    = $i;
        }
    }

    if($min_distance == -1) {
        return;
    }

    return $best_point;
}
=cut

sub get_most_remote_node {
    my $self = shift;
    my $cluster_index = shift;

    my $c = $self->centroids->[$cluster_index];

    my $max_d = 0;
    my $max_id = -1;
    foreach my $assoc (@{$self->associations}) {
        next unless $assoc->[1] == $cluster_index;
        my $p = $self->points->[$assoc->[0]];
        my $d = $self->get_distance($c,$p);
        if($max_d < $d) {
            $max_d = $d;
            $max_id = $p;
        }
    }

    return unless $max_id != -1;
    return $max_id;
}

sub get_cluster_radius {
    my $self = shift;
    my $cluster_index = shift;

    my $c = $self->centroids->[$cluster_index];

    my $max_d = 0;
    foreach my $assoc (@{$self->associations}) {
        next unless $assoc->[1] == $cluster_index;
        my $p = $self->points->[$assoc->[0]];
        my $d = $self->get_distance($c,$p);
        if($max_d < $d) {
            $max_d = $d;
        }
    }

    return 1 unless $max_d != 0;
    return $max_d;
}

#takes a node index and a centroid index.
#mass * location * density.
sub move_centroid {
    my($self, $centroid_index, $point_radius) = @_;

    my @points;
    foreach my $assoc (@{$self->associations}) {
        push @points, $assoc->[0] if($assoc->[1] == $centroid_index);
    }

    my $x = 0;
    my $y = 0;
    my $m = 0;
    foreach my $point_index (@points) {
        my $mass = $self->get_point_density($point_index, $point_radius);
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
=pod
    my($self, $centroid_index, $radius) = @_;

    my $total_mass = 0;
    my @sums = (0,0);
    for(my $i = 0 ; $i < scalar @{$self->points} ; $i++) {
        if($self->get_assoc_for_pt($i) == $centroid_index) {
            $sums[0] += $self->get_point_density($i , $centroid_index, $radius) *
            $self->points->[$i]->[0];
            
            $sums[1] += $self->get_point_density($i , $centroid_index, $radius) *
            $self->points->[$i]->[0];  

            $total_mass += $self->get_point_density($i , $centroid_index,$radius);
        }
    }

    if($total_mass != 0) {
        $sums[0] = $sums[0] / $total_mass;
        $sums[1] = $sums[1] / $total_mass;
    }

    $self->centroids->[$centroid_index] = \@sums; 
}
=cut

=pod

    my($self , $centroid_index) = @_;
    cluck 'wtf' unless $self->centroids->[$centroid_index];

    print "\n\n$centroid_index info.\n";
    my $radius = $self->get_cluster_radius($centroid_index);
    print "Using radius $radius\n";

    #this code will find the sum of the coordinates and their densities,
    #comprising the numerator of My = sum[p(x) * m sub x] / m
    my @center = (0,0);
    for(my $i = 0 ; $i < scalar @{$self->points} ; $i++) {
        if($self->get_assoc_for_pt($i) == $centroid_index) {
            #p.sure this is the actually right density :)
            my $density  = $self->get_node_density($i,$centroid_index,$radius);

            print "using density at $i as $density\n";

            $center[0] += $density * $self->points->[$i]->[0];
            $center[1] += $density * $self->points->[$i]->[1];
        }
    }

    print "cetner: ";
    print Dumper \@center;


    #this code will divide by the total mass of the the circle
    #note: get cluster density isn't dependent on anything except the number of
    #nodes in the cluster remaining constant, which they do here.
    #
    print "starting point density for cluster $centroid_index is: " .
    $self->get_cluster_density($centroid_index) . "\n";

    my $total_mass = $self->_get_total_mass($radius ,
        $self->get_cluster_density($centroid_index)); 

    print "total mass for $centroid_index is $total_mass\n";

    $center[0] = $center[0] / $total_mass;
    $center[1] = $center[1] / $total_mass;

    $self->centroids->[$centroid_index] = \@center; 

    print Dumper \@center;


    print "\n\n";
}
=cut

#density of a node = point's density + centroid's lamina density at that point.
#centroid's lam den = (-centroid density / radius) * distance from pt to
#centroid
sub get_node_density {
    my($self , $point_index , $centroid_index , $radius) = @_;

    my $p_d = $self->get_point_density($point_index,$radius);
    my $c_d = $self->get_cluster_density($centroid_index);
    my $distance = $self->get_distance($self->points->[$point_index] ,
    $self->centroids->[$centroid_index]);

    return (((-$c_d/$radius) * ($distance)) + $c_d);
}

sub move_centroids {
    my($self, $radius) = @_;
    for(my $i=0;$i<scalar @{$self->centroids};$i++) {
        $self->move_centroid($i, $radius);
    }
}

#This code finds the total mass of a cluster.
#using simpson's rule w/ n=20 to estimate integral b/c it's getting late.
sub _get_total_mass {
    my($self , $radius , $cluster_density) = (undef , 1 , 2);

    #this is the function to integrate. just see the paper...
    #note: this is broken b/c mue depends on x and y, so (-c/r)*x != true mue.
    my $f = sub {
        my($r , $c , $x) = @_; #radius , cluster_density, x
        return sqrt($r*$r - ($x*$x - 2*$r*$x + $r*$r)) * ((-$c/$r) * $x + $c);
    };


    my $term_sum = 0;
    for(my $i = 0 ; $i <= 50 ; $i++) {
        my $coeff = 1;
        if($i == 0 || $i == 50) {
            $coeff = 1;
        }
        elsif($i % 2 == 0) {
            $coeff = 2;
        }
        else {
            $coeff = 4;
        }
        $term_sum += $coeff * $f->($radius , $cluster_density , $i/50);
    }

    $term_sum = $term_sum * 2; #add bottom half of circle back in.

    return (((2*$radius)/(3*50)) * $term_sum);
}


#this will run a single pass of the clustering.
sub force_pass {
    my($self, $radius) = @_;

    my $f = `ls /home/nfulton/cs/geometry/imgz/ | grep force | wc -l`;
    chop $f;
    $self->get_r_output("force-$f");

    #print "Moving Clusters\n";

    for(my $i=0;$i<scalar @{$self->centroids};$i++) {
        my $closest = $self->find_closest_point($i);
        if(defined $closest and $closest != -1) {
            print "doing stuff\n";
            $self->set_association($closest, $i);
        }
        else {
            print "$closest is null.\n" if defined $closest;
        }
    }

    my @old_centroids = @{$self->centroids};
    $self->move_centroids($radius);
}

sub run_pass {
    my($self , $radius , $disable_images) = @_;
        
    my $count = 0;

    my $changed = 1;
    while($changed or scalar @{$self->points} > scalar @{$self->associations}) {
        $self->get_r_output($count) unless defined $disable_images;
        $count++;

        #print "Moving Clusters\n";

        $changed = 0;
        for(my $i=0;$i<scalar @{$self->centroids};$i++) {
            my $closest = $self->find_closest_point($i);
            if(defined $closest and $closest != -1) {
                $self->set_association($closest, $i);
                $changed = 1;
            }
        }


        my @old_centroids = @{$self->centroids};
        $self->move_centroids($radius);


#        for(my $j = 0 ; $j < scalar @old_centroids ; $j++) {
#            my($x1,$y1) = @{$old_centroids[$j]};
#            my($x2,$y2) = @{$self->centroids->[$j]};
#
#            ($x1, $x2, $y1,$y2) = ( nearest(.001, $x1) , nearest(.001, $x2),
#            nearest(.001, $y1),
#            nearest(.001, $y2));
#
#            if($x1 != $x2 || $y1 != $y2) {
#                $changed = 1;
#            }
#
#        }
    }

    #print "Done Clustering...\n";
}

=pod
    #there's nothing adaptive here... it's actually pretty crappy.
    #just run through all of the points and set their association.
    for(my $runs = 0 ; $runs < scalar @{$self->points} / scalar @{$self->centroids} ; $runs++) {
        for(my $i = 0 ; $i < scalar @{$self->centroids} ; $i++) {
            my $closest = $self->find_closest_point($i);
            if(defined $closest and $closest != -1) {
                print "$i choosing point $closest\n";
                $self->set_association($closest , $i);
            }
            else {
                print "$i choosing no point.\n";
            }

            #$self->print_results();
        }

        $self->move_centroids($radius);
    }
}
=cut


sub get_r_output {
    print "Generating R output...\n";
    my $self = shift;
    my $file_name = shift;
    my @colors = ('green','blue','red','purple','darkgoldenrod','aquamarine');


    my $biggest_x=0;
    my $biggest_y=0;
    foreach my $p (@{$self->points}) {
        $biggest_x = $p->[0] if $p->[0] > $biggest_x;
        $biggest_y = $p->[1] if $p->[1] > $biggest_y;
    }
    $biggest_x += 6;
    $biggest_y += 6;
    my $plot = "png(\"/home/nfulton/cs/geometry/imgz/$file_name.png\", bg=\"white\", width=700, height=500);plot(NULL,xlim=c(0,$biggest_x),ylim=c(0,$biggest_y),xlab=\"x\",ylab=\"y\",main=\"Clustering Results\");";


    my $clusters = scalar @{$self->centroids};


    my $cmd = "";

    my @x_cluster = ();
    my @y_cluster = ();
    #for(my $i=0;$i<$clusters;$i++) { $x_cluster[$i] = "" ; $y_cluster[$i] = "";}

    foreach my $assoc (@{$self->associations}) {
        $x_cluster[$assoc->[1]] .= $self->points->[$assoc->[0]]->[0] . ",";
        $y_cluster[$assoc->[1]] .= $self->points->[$assoc->[0]]->[1] . ",";
    }

    for(my $i=0;$i<scalar @x_cluster;$i++) {
        my $x = $x_cluster[$i];
        my $y = $y_cluster[$i];
        next unless $x and $y;
        chop $x; chop $y;

        $x = "cx$i <- c(" . $x . ");";
        $y = "cy$i <- c(" . $y . ");";

        $cmd .= $x . $y . " points(cx$i, cy$i,col=\"".$colors[$i]."\");";
    }


    my $ret = "";
    for(my $i=0;$i<scalar @{$self->centroids};$i++) {
        my $x = $self->centroids->[$i]->[0];
        my $y = $self->centroids->[$i]->[1];
        my $label = substr $colors[$i] , 0, 1;

        $ret .= "text($x,$y,labels=\"$label\");";
    }
=pod
    my $ret = "c_x <- c(";
    foreach my $x (@{$self->centroids}) {
        $ret .= $x->[0] . ",";
    }
    chop $ret;
    $ret .= ") ; c_y <- c(";
    foreach my $y (@{$self->centroids}) {
        $ret .= $y->[1] . ",";
    }
    chop $ret;

    $ret .= ");points(c_x,c_y,col=\"".$colors[0]."\");\n";
=cut

    my $empty = "";
    for(my $pi = 0; $pi < scalar @{$self->points} ; $pi++) {
        if($self->get_assoc_for_pt($pi) == -1) {
            my $x = $self->points->[$pi]->[0];
            my $y = $self->points->[$pi]->[1];
            $empty .= "points($x,$y,col=\"black\");";
        }
    }



    my $R = Statistics::R->new() ;
    $R->startR();
    $R->send($plot . $cmd .  $empty . $ret . "dev.off();");
    $R->stopR();
}


#prints output.
sub print_results {
    my($self) = @_;
    for(my $i = 0 ; $i < scalar @{$self->centroids} ; $i++) {
        print "Centroid Number " . $i . "\n";
        print "LOCATION OF CENTROID: " . $self->centroids->[$i]->[0] . "," .
        $self->centroids->[$i]->[0] . "\n";
        print "And the points: \n";


        for(my $l = 0; $l < scalar @{$self->points} ; $l++) {
            my @point = @{$self->points->[$l]};
            if($self->get_assoc_for_pt($l) == $i) {
                print $point[0] . "," . $point[1] . "\n";
            }
        }

    }

}


sub davies_bouldin {
    my $self = shift;
    
    my $max = -1;

    my $n = scalar @{$self->centroids};
    for(my $i = 0; $i < $n ; $i++) {
        for(my $j = 0; $j < $n ; $j++) {
            next if $i == $j;
            my $d = $self->get_distance($self->centroids->[$i],
                $self->centroids->[$j]);
            my $f = ($self->_gad($i) + $self->_gad($j)) / $d;

            if($f > $max || $max == -1) {
                $max = $f;
            }
        }
    }

    return $max;
}
#get avg distance
sub _gad {
    my($self, $centroid_index) = @_;
    my $c= $self->centroids->[$centroid_index];

    my $d = 0;
    my $t = 0;
    foreach $a (@{$self->associations}) {
        next if $a->[1] != $centroid_index;
        my $p = $self->points->[$a->[0]] ;
        $d += $self->get_distance($c,$p);
        $t++;
    }

    return $d/$t;
}

1;


for(my $dbi_count = 0 ; $dbi_count < 20 ; $dbi_count++) {
    my $c = new Cluster;

    my $rand = 200;
    for(my $i=0;$i<$rand;$i++) {
        my $x = int(rand(20));
        my $y = int(rand(20));

        $c->add_point($x,$y);
    }

#$c->add_point(1,1);
#$c->add_point(1,2);
#$c->add_point(1,5);
#$c->add_point(2,2);
#$c->add_point(3,3);
#$c->add_point(7,1);
#$c->add_point(7,3);
#$c->add_point(8,2);
#$c->add_point(8,5);
#$c->add_point(9,1);
#$c->add_point(10,3);
#$c->add_point(11,2);
#$c->add_point(11,8);


    $c->add_centroid(0,0);
    $c->add_centroid(0,0);
    $c->add_centroid(0,0);
    $c->add_centroid(0,0);
    $c->add_centroid(0,0);

    $c->run_pass(2,1);
#$c->get_r_output("final");

    print "Davies-Bouldin Index: " . $c->davies_bouldin . "\n";
}
