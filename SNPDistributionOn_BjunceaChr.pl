#!usr/bin/perl -w
use SVG;

my($Width, $Height) = (10000, 6000);
my $svg = SVG -> new(width => $Width, height => $Height);
my %chr = ();
open IN0, $ARGV[0]; ##- chromosome sizes
while(<IN0>){
  chomp;
  my ($chr, $length) = (split(/\t/, $_))[0,1];
  if(not exists $chr{$chr}){
     $chr{$chr} = $length;
  }
  else{
     die "Error code 1.\n";
  }
}
close IN0;


my $WaneRatio = 7000/70000000; #### one base pair occupied pixels

my $x_origin = $Width - 8500; #origin position of x axis
my $y_origin = $Height - 5500; #origin position of y axis

$svg -> line(x1 => $x_origin, y1 => $y_origin - 50, x2 => $x_origin+7500, y2 => $y_origin - 50, stroke => 'black', 'stroke-width' =>20);

my $Size = 7000/10;
my $Value_y = 7;

for(my $k = 0; $k <= 10; $k++)
{
        $svg -> line(y1 => $y_origin - 50, x1 => $x_origin + $Size*$k, y2 => $y_origin - 100, x2 => $x_origin + $Size*$k, stroke => 'black','stroke-width' => 10 );
}

####Draw number of scales on y axis
for(my $k = 0; $k <= 10; $k++)
{
	my $Number = $Value_y * $k;
        $svg -> text(x => $x_origin + $Size*$k - 170, y => $y_origin + 120, 'font-size' => 150, stroke => 'black', -cdata => "$Number Mb");
}


####Draw chromosome
my $Track_size = 280;
my $Chr1_Position = $y_origin + 350;
my $Chr1_Position_Y = $y_origin + 320;
my @allchr = sort keys %chr;
for (my $k = 0; $k <= $#allchr; $k++){
     $svg -> text(x => $x_origin - 250, y => $Chr1_Position + $Track_size*$k + 100, 'font-size' => 100, stroke => 'black', -cdata => "$allchr[$k]");  
     $svg -> rect(x => $x_origin, y => $Chr1_Position + $Track_size*$k,  width => ($chr{$allchr[$k]})*$WaneRatio, height => 150, fill => "lightblue");
      

     ##- read gene file
     open IN1, $ARGV[1];
     while(<IN1>){
        chomp;
        my @temp = split(/\t/, $_);
        if($allchr[$k] eq $temp[1]){
         #      $svg -> rect(x => $x_origin+$temp[2]*$WaneRatio, y => $Chr1_Position + $Track_size*$k,  width => ($temp[3] - $temp[2])*$WaneRatio, height => 150, fill => "red");
                $svg -> rect(x => $x_origin+$temp[2]*$WaneRatio-5, y => $Chr1_Position + $Track_size*$k,  width => 11, height => 150, fill => "red");
        }
     }
     close IN1;

     ##- read snp info
     open IN2, $ARGV[2];
     while(<IN2>){
       chomp;
       my @temp = split(/\t/, $_);
       if($allchr[$k] eq $temp[0]){
       #   $svg -> line(x1 => $x_origin + $temp[1]*$WaneRatio, y1 => $Chr1_Position + $Track_size *$k - 20, 
       #                x2 => $x_origin + $temp[1]*$WaneRatio, y2 => $Chr1_Position + $Track_size *$k, stroke=>'blue',"stroke-width"=>20);
           my ($x1, $y1) = ($x_origin + $temp[1]*$WaneRatio, $Chr1_Position + $Track_size*$k);
           my ($x2, $y2) = ($x_origin + $temp[1]*$WaneRatio-20, $Chr1_Position + $Track_size*$k-30);
           my ($x3, $y3) = ($x_origin + $temp[1]*$WaneRatio+20, $Chr1_Position + $Track_size*$k-30);
           $svg -> polygon(points=>"$x1,$y1 $x2,$y2 $x3,$y3", style=>"fill:blue;stroke:blue;stroke-width:1");
       }

     }
     close IN2;

     ##- read indel
     open IN3, $ARGV[3];
     while(<IN3>){
       chomp;
       my @temp = split(/\t/, $_);
       if($allchr[$k] eq $temp[0]){
       #  $svg -> line(x1 => $x_origin + $temp[1]*$WaneRatio, y1 => $Chr1_Position + $Track_size *$k + 150,
       #               x2 => $x_origin + $temp[1]*$WaneRatio, y2 => $Chr1_Position + $Track_size *$k + 170, stroke=>'lime',"stroke-width"=>20);
          my ($x1, $y1) = ($x_origin + $temp[1]*$WaneRatio, $Chr1_Position + $Track_size *$k + 150);
          my ($x2, $y2) = ($x_origin + $temp[1]*$WaneRatio-20, $Chr1_Position + $Track_size *$k + 150+30);
          my ($x3, $y3) = ($x_origin + $temp[1]*$WaneRatio+20, $Chr1_Position + $Track_size *$k + 150+30);
          $svg -> polygon(points=>"$x1,$y1 $x2,$y2 $x3,$y3", style=>"fill:lime;stroke:lime;stroke-width:1");
       }
     }
     close IN3;
     
     ##- read stop snp
     open IN4, $ARGV[4];
     while(<IN4>){
        chomp;
        my @temp = split(/\t/, $_);
        if($allchr[$k] eq $temp[0]){
           my ($x1, $y1) = ($x_origin + $temp[1]*$WaneRatio, $Chr1_Position + $Track_size*$k);
           my ($x2, $y2) = ($x_origin + $temp[1]*$WaneRatio-20, $Chr1_Position + $Track_size*$k-30);
           my ($x3, $y3) = ($x_origin + $temp[1]*$WaneRatio+20, $Chr1_Position + $Track_size*$k-30);
           $svg -> polygon(points=>"$x1,$y1 $x2,$y2 $x3,$y3", style=>"fill:red;stroke:red;stroke-width:1");
        }    
     }
     close IN4;

}

open OUT,">$ARGV[5]";
print OUT $svg->xmlify;
system("/home/caix/bin/distributing_svg_4.74/svg2xxx_release/svg2xxx  $ARGV[5] -t pdf");

