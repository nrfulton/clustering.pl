=pod

Cluster is a role that a Cluster::Algorithm implements.

Cluster provides the following data:
    centroids    - each value is a refernce to [x][y]
    points       - each value is a reference to [x][y]
    assoications - each value is a [point_index][centroid_index]
    make_images  - if true images will be generated in do_cluster(); false by
                   default.

Cluster provides access to the following methods; they can be overloaded:
    make_r_img(image_location)  - creates a PNG file at image_location.
                                  must be an absolute path!!!
    get_distance([x][y],[x][y]) - returns the distance between the two points.
    do_cluster() - does do_iteration until do_iteration returns two successive
                   true values.


Cluster requires the following methods to be implemented by comsumers:
    do_iteration - does a single iteration of a cluster algorithm
                   should return a true value if the cluster has converged.
    

=cut

use strict;
use warnings;

package Cluster; 
use Moose::Role;
use Data::Dumper;
use  Carp qw(cluck);
use Statistics::R;
use Math::Round;

requires 'do_iteration';

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
    sub add_point {
        my $self = shift;
        push @{$self->points} , \@_;
    }

#A list of the coordinates of centroids.
has 'centroids' => (
    is  => 'rw',
    isa => 'ArrayRef',
    default => sub {
        my @array;
        return \@array;
    },
);
    sub add_centroid {
        my $self = shift;
        push @{$self->centroids} , \@_;
    }

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

has 'make_images' => (
    is => 'rw',
    default => sub { return 0; }
);

has 'do_cluster_count' => (
    is => 'rw',
    default => sub { return 0; }
);

has 'output_directory' => (
    is => 'rw',
    default => sub { return "/tmp"; }
);

sub get_distance {
    my($self,$c1,$c2)=@_;
    cluck('err') unless ref $c1 eq "ARRAY" and ref $c2 eq "ARRAY";
    return sqrt(abs(($c2->[1] - $c1->[1]) + ($c2->[0] - $c1->[0])));
}



#The following methods control R via the Perl/R bridge, allowing access to 
sub make_r_image {
    print "Generating R input...\n";

    my $self = shift;
    my $step_count = shift;
    my @colors = ('green','blue','red','purple','darkgoldenrod','aquamarine');

    my $R = Statistics::R->new() ;
    $R->startR();
    $R->send($self->_make_device_cmd($step_count).
             $self->_make_plot_cmd . 
             $self->_make_points_cmd . 
             $self->_make_centroids_cmd(\@colors) . 
             $self->_make_nodes_cmd(\@colors) . 
             "dev.off();");
    $R->stopR();
    return 1;
}
    sub _make_plot_cmd {
        my $self = shift;

        #The xmax and ymax values for the plot must be set when the plot
        #cmd is called.
        my $biggest_x=0;
        my $biggest_y=0;
        foreach my $p (@{$self->points}) {
            $biggest_x = $p->[0] if $p->[0] > $biggest_x;
            $biggest_y = $p->[1] if $p->[1] > $biggest_y;
        }
        foreach my $p (@{$self->centroids}) {
            $biggest_x = $p->[0] if $p->[0] > $biggest_x;
            $biggest_y = $p->[1] if $p->[1] > $biggest_y;
        }

        ($biggest_x,$biggest_y) = ($biggest_x + 1 , $biggest_y + 1);
        return "plot(NULL,xlim=c(0,".$biggest_x."), ".
               "ylim=c(0,".$biggest_y."), xlab=\"x\", ylab=\"y\", ".
               "main=\"Clustering Results\");";
    }
 
    sub _make_device_cmd {
        my($self, $step_count) = @_;
        my $save_dir = $self->_get_full_directory;
        `mkdir -p $save_dir` unless -e $save_dir;
        print "Saving to $save_dir/$step_count.png";
        return "png(\"$save_dir/$step_count.png\", ".
               "bg=\"white\", width=700, heigh=500);";
    }

    sub _get_full_directory {
        my($self) = @_;
        my(undef,undef,undef,undef,undef,undef,$year,undef,$year_day,undef) =
        localtime(time);
        return $self->output_directory . "/clustering_output/" . 
        $year ."-".$year_day .
                scalar @{$self->points} ."-". scalar @{$self->centroids} ."-". $$
                ."/".$self->do_cluster_count;
    }

    sub _get_clusters {
        my($self) = @_;
        my @x_cluster = ();
        my @y_cluster = ();
    
        foreach my $assoc (@{$self->associations}) {
            $x_cluster[$assoc->[1]] .= $self->points->[$assoc->[0]]->[0] . ",";
            $y_cluster[$assoc->[1]] .= $self->points->[$assoc->[0]]->[1] . ",";
        }

        return(\@x_cluster , \@y_cluster);
    }

    sub _make_nodes_cmd { #x_cluster = list of cluster's x coords
        my($self, $colors) = @_;
        my($x_cluster, $y_cluster) = $self->_get_clusters;

        my $ret_var = "";
        for(my $i=0;$i<scalar @{$x_cluster};$i++) {
            my $x = $x_cluster->[$i];
            my $y = $y_cluster->[$i];
            next unless $x and $y;
            chop $x; chop $y;

            $x = "cx$i <- c(" . $x . ");";
            $y = "cy$i <- c(" . $y . ");";

            $ret_var .= $x . $y . 
                        " points(cx$i, cy$i,col=\"".$colors->[$i]."\");";
        }

        return $ret_var;
    }
    sub _make_points_cmd {
        my($self) = @_;

        my $empty = "";
        for(my $pi = 0; $pi < scalar @{$self->points} ; $pi++) {
            if($self->get_assoc_for_pt($pi) == -1) {
                my $x = $self->points->[$pi]->[0];
                my $y = $self->points->[$pi]->[1];
                $empty .= "points($x,$y,col=\"black\");";
            }
        }

        return $empty;
    }
    sub _make_centroids_cmd {
        my $self = shift;
        my $colors = shift;
        my $ret = "";
        for(my $i=0;$i<scalar @{$self->centroids};$i++) {
            my $x = $self->centroids->[$i]->[0];
            my $y = $self->centroids->[$i]->[1];
            my $label = substr $colors->[$i] , 0, 1;

            $ret .= "text($x,$y,labels=\"$label\");";
        }
        return $ret;
    }


sub do_cluster {
    my $self = shift;

    my $truth_counter = 0;

    my $pid = $$;

    my $step_count = 0;
    while($truth_counter != 2) {
        if($self->make_images) {
            $self->make_r_image($step_count);
        }
        $truth_counter += 1 if $self->do_iteration(@_); #pass the stack along.

        $step_count+=1;
    }

    $self->do_cluster_count($self->do_cluster_count + 1);
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

sub get_assoc_for_pt {
    my $self = shift;
    my $pt = shift;

    map {
        return $_->[1] if $_->[0] == $pt;
    } @{$self->associations};

    return -1;
}



1;

