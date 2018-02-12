use lib 'lib';
use Cluster::Inn;



for(my $dbi_count = 0 ; $dbi_count < 1 ; $dbi_count++) {
    my $c = new Cluster::Inn;

    my $pts_num = 20; #the number of points you want in the data set to cluster.
    for(my $i=0;$i<$pts_num;$i++) {
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


    $c->make_images(1);

    $c->add_centroid(0,0);
    $c->add_centroid(0,0);
    $c->add_centroid(0,0);
    $c->add_centroid(0,0);
    $c->add_centroid(0,0);

    $c->do_cluster(2);
#$c->get_r_output("final");

    print "Davies-Bouldin Index: " . $c->davies_bouldin . "\n";
}

