use strict;
use Bio::SearchIO; 

my $fileIn = $ARGV[0];
my $saveCnt = $ARGV[1];

my $in = new Bio::SearchIO(-format => 'blast', 
                           -file   => $fileIn);
                           
open (FILEOUT, ">$saveCnt");       

                 
                           
while( my $result = $in->next_result ) 
{
  ## $result is a Bio::Search::Result::ResultI compliant object
     my $query_id  =  $result->query_name;
     my $query_len =  $result->query_length;
    
    #
    
    my $group = 0;
    while( my $hit = $result->next_hit ) {
    	
    ## $hit is a Bio::Search::Hit::HitI compliant object
         my $hit_id  =   $hit->name;
         my $hit_len =   $hit->length;
         my $desc    = $hit->description;
   
         $group++;

         while( my $hsp = $hit->next_hsp ) {
              ## $hsp is a Bio::Search::HSP::HSPI compliant object
              #my $hsp_len = $hsp->num_identical;
               my $score =  $hsp->score ;
               my $evalue=  $hsp->evalue;
               #my @query_range = $hsp->range('query');
               #my @hit_range   = $hsp->range('hit');
               
               my $start_q    = $hsp->	start('query') ;
               my $end_q      = $hsp->	end('query') ;
               my $start_h    = $hsp->	start('hit') ;
               my $end_h      = $hsp->	end('hit');
               
               my $strand     = $hsp->query->strand;
               
               
               

               my $each_group = $group;
               
               my @output = ($query_id,$query_len,$start_q, $end_q  ,$each_group,$hit_id,$hit_len ,$start_h, $end_h,$score,$evalue,$strand,$desc);
               
               print FILEOUT join("\t",@output),"\n";
             
         }

   }  
}
