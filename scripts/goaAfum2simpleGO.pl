use strict;
use warnings;

# process GOA file 
# eg ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/22118.N_fumigata_ATCC_MYA-4609.goa
# for use in bioconductor R analysis of GO enrichment (need a table of gene to GO term)

while(<>){
    next if /^\!/;
    chomp;
    my @row = split(/\t/,$_);
#    my $i = 0;
#    for my $n ( @row ) {
#	print "$i=$n\n";
#	$i++;
#    }
    my $go = $row[4];
    my $evid = $row[6];
    my $genename = $row[2];
    my $full_ids = $row[10];
    my $systematicname = $genename;
    $systematicname =~ s/^(AFUA_(\S+))/$1\tAfu\L$2/;
    my $alias_name = $genename;
    for my $id ( split(/\|/,$full_ids) ) {
	$id =~ s/\/(\d+)$//;
	if( $id =~ /AFUA_(\S+)/ ) {
	    $alias_name = $id;
	    $systematicname = sprintf("Afu%s",lc($1));
	    last;
	}
    }
    print join("\t", $systematicname,$evid,$go,$alias_name,$genename),"\n";
}
