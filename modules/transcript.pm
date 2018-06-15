package transcript;
use strict;



#nom, second_name, chr, brin, nbre_exons, nom_prot, short_prot, taille_prot, uniprot_id, acc_version, gi_NM, acc_g, gi_NG, acc_p, gi_NP, ENST, ENSP,  translation_start_site, mutalyzer_version, main, USHER, RP, DFN, "MiSeq-28", "MiSeq-112", "MiSeq-121", "MiSeq-132", "MiSeq-3", "MiniSeq-3",  "MiniSeq-121", "MiniSeq-132"
#need new attributes EXTENDED_DFN et "NextSeq-ClinicalExome"

sub new { #constructeur 
        my ($class, $nm) = @_;
        my $self;
        if (!defined($nm)) {die "you must provide a name (NM)!!!"}
        else {
		$self = {
                        'nm' => $nm,
			'second_name' => 'NULL',
			'ng' => 'NG_000000.0',
			'np' => 'NP_000000.0',
			'nm_version' => '1',
			'main' => 'f',
			'mutal' => '_v001',
                        'USHER' => 'f',
                        'RP' => 'f',
                        'DFN' => 'f',
                        'EXTENDED_DFN' => 't',
			'"MiSeq-28"' => 'f',
			'"MiSeq-112"' => 'f',
			'"MiSeq-121"' => 'f',
			'"MiSeq-132"' => 'f',
			'"MiSeq-3"' => 'f',
			'"MiniSeq-3"' => 'f',
			'"MiniSeq-121"' => 'f',
			'"MiniSeq-132"' => 'f',
			'"MiniSeq-152"' => 't',
			'"NextSeq-ClinicalExome"' => 't',
		};
        }
	bless ($self, $class);
	return $self;
}

sub setNM {
        my ($self, $nm) = @_;
        $self->{nm} = $nm;
}
sub getNM {
        my $self = shift;
        return $self->{nm};
}
sub setGeneName {my ($self, $gene_name) = @_;$self->{gene_name} = $gene_name}
sub getGeneName {my $self = shift;return $self->{gene_name}}
sub setSecondName {my ($self, $second_name) = @_;$self->{second_name} = $second_name}
sub getSecondName {my $self = shift;return $self->{second_name}}
sub setChr {my ($self, $chr) = @_;$self->{chr} = $chr}
sub getChr {my $self = shift;return $self->{chr}}
sub setStrand {my ($self, $strand) = @_;$self->{strand} = $strand}
sub getStrand {my $self = shift;return $self->{strand}}
sub setNbExons {my ($self, $nb_exons) = @_;$self->{nb_exons} = $nb_exons}
sub getNbExons {my $self = shift;return $self->{nb_exons}}
sub setProtName {my ($self, $prot_name) = @_;$self->{prot_name} = $prot_name}
sub getProtName {my $self = shift;return $self->{prot_name}}
sub setShortProt {my ($self, $short_prot) = @_;$self->{short_prot} = $short_prot}
sub getShortProt {my $self = shift;return $self->{short_prot}}
sub setProtSize {my ($self, $prot_size) = @_;$self->{prot_size} = $prot_size}
sub getProtSize {my $self = shift;return $self->{prot_size}}
sub setUniprot {my ($self, $uniprot) = @_;$self->{uniprot} = $uniprot}
sub getUniprot {my $self = shift;return $self->{uniprot}}
sub setNMVersion {my ($self, $nm_version) = @_;$self->{nm_version} = $nm_version}
sub getNMVersion {my $self = shift;return $self->{nm_version}}
sub setgiNM {my ($self, $ginm) = @_;$self->{ginm} = $ginm}
sub getgiNM {my $self = shift;return $self->{ginm}}
sub setNG {my ($self, $ng) = @_;$self->{ng} = $ng}
sub getNG {my $self = shift;return $self->{ng}}
sub setgiNG {my ($self, $ging) = @_;$self->{ging} = $ging}
sub getgiNG {my $self = shift;return $self->{ging}}
sub setNP {my ($self, $np) = @_;$self->{np} = $np}
sub getNP {my $self = shift;return $self->{np}}
sub setgiNP {my ($self, $ginp) = @_;$self->{ginp} = $ginp}
sub getgiNP {my $self = shift;return $self->{ginp}}
sub setENST {my ($self, $enst) = @_;$self->{enst} = $enst}
sub getENST {my $self = shift;return $self->{enst}}
sub setENSP {my ($self, $ensp) = @_;$self->{ensp} = $ensp}
sub getENSP {my $self = shift;return $self->{ensp}}
sub setTss {my ($self, $tss) = @_;$self->{tss} = $tss}
sub getTss {my $self = shift;return $self->{tss}}
sub setMutalyzerVer {my ($self, $mutal) = @_;$self->{mutal} = $mutal}
sub getMutalyzerVer {my $self = shift;return $self->{mutal}}
sub setMain {my ($self, $main) = @_;$self->{main} = $main}
sub getMain {my $self = shift;return $self->{main}}
sub setUsh {my ($self, $ush) = @_;$self->{USHER} = $ush}
sub getUsh {my $self = shift;return $self->{USHER}}
sub setRP {my ($self, $rp) = @_;$self->{RP} = $rp}
sub getRP {my $self = shift;return $self->{RP}}
sub setDfn {my ($self, $dfn) = @_;$self->{DFN} = $dfn}
sub getDfn {my $self = shift;return $self->{DFN}}
sub setExtDfn {my ($self, $ext_dfn) = @_;$self->{EXTENDED_DFN} = $ext_dfn}
sub getExtDfn {my $self = shift;return $self->{EXTENDED_DFN}}
sub setMS28 {my ($self,$ms28 ) = @_;$self->{'"MiSeq-28"'} = $ms28}
sub getMS28 {my $self = shift;return $self->{'"MiSeq-28"'}}
sub setMS112 {my ($self, $ms112) = @_;$self->{'"MiSeq-112"'} = $ms112}
sub getMS112 {my $self = shift;return $self->{'"MiSeq-112"'}}
sub setMS121 {my ($self, $ms121) = @_;$self->{'"MiSeq-121"'} = $ms121}
sub getMS121 {my $self = shift;return $self->{'"MiSeq-121"'}}
sub setMS132 {my ($self, $ms132) = @_;$self->{'"MiSeq-132"'} = $ms132}
sub getMS132 {my $self = shift;return $self->{'"MiSeq-132"'}}
sub setMS3 {my ($self, $ms3) = @_;$self->{'"MiSeq-3"'} = $ms3}
sub getMS3 {my $self = shift;return $self->{'"MiSeq-3"'}}
sub setMnS3 {my ($self, $mns3) = @_;$self->{'"MiniSeq-3"'} = $mns3}
sub getMnS3 {my $self = shift;return $self->{'"MiniSeq-3"'}}
sub setMnS121 {my ($self, $mns121) = @_;$self->{'"MiniSeq-121"'} = $mns121}
sub getMnS121 {my $self = shift;return $self->{'"MiniSeq-121"'}}
sub setMnS132 {my ($self, $mns132) = @_;$self->{'"MiniSeq-132"'} = $mns132}
sub getMnS132 {my $self = shift;return $self->{'"MiniSeq-132"'}}
sub getMnS152 {my $self = shift;return $self->{'"MiniSeq-152"'}}
sub setNxTSqCE {my ($self, $nxtsq_ce) = @_;$self->{'"NextSeq-ClinicalExome"'} = $nxtsq_ce}
sub getNxTSqCE {my $self = shift;return $self->{'"NextSeq-ClinicalExome"'}}
sub setHGNC {my ($self, $hgnc) = @_;$self->{hgnc} = $hgnc}
sub getHGNC {my $self = shift;return $self->{hgnc}}


sub toPrint {
	my $self = shift;
	my $txt =  "#Gene\tSecondName\tNM\tENST\tENSP\tHGNC\tUNIPROT\tNMVersion\tNG\tNP\tMainIsoform\tChr\tStrand\t#Exons\tTssPos\tpName\tpShort\tpSize\n";
	$txt .= $self->getGeneName()."\t".$self->getSecondName()."\t".$self->getNM()."\t".$self->getENST()."\t".$self->getENSP()."\t".$self->getHGNC()."\t".$self->getUniprot()."\t".$self->getNMVersion()."\t".$self->getNG()."\t".$self->getNP()."\t".$self->getMain()."\t".$self->getChr()."\t".$self->getStrand()."\t".$self->getNbExons()."\t".$self->getTss()."\t".$self->getProtName()."\t".$self->getShortProt()."\t".$self->getProtSize()."\n";
	return $txt;
}

sub toSQL {
	my $self = shift;
	my $chr = $self->getChr();
	$chr =~ s/chr//o;
	my $sql = "INSERT INTO gene (nom, second_name, chr, brin, nbre_exons, nom_prot, short_prot, taille_prot, uniprot_id, acc_version, gi_nm, acc_g, gi_ng, acc_p, gi_np, enst, ensp,  translation_start_site, mutalyzer_version, main, usher, rp, dfn, extended_ns, \"MiSeq-28\", \"MiSeq-112\", \"MiSeq-121\", \"MiSeq-132\", \"MiSeq-3\", \"MiniSeq-3\", \"MiniSeq-121\", \"MiniSeq-132\", \"NextSeq-ClinicalExome\", \"MiniSeq-152\") VALUES ( '{\"".$self->getGeneName()."\",\"".$self->getNM()."\"}','".$self->getSecondName()."','".$chr."','".$self->getStrand()."','".$self->getNbExons()."','".$self->getProtName()."','".$self->getShortProt()."','".$self->getProtSize()."','".$self->getUniprot()."','".$self->getNMVersion()."',NULL,'".$self->getNG()."',NULL,'".$self->getNP()."',NULL,'".$self->getENST()."','".$self->getENSP()."','".$self->getTss()."','".$self->getMutalyzerVer()."','".$self->getMain()."','".$self->getUsh()."','".$self->getRP()."','".$self->getDfn()."','".$self->getExtDfn()."','".$self->getMS28()."','".$self->getMS112()."','".$self->getMS121()."','".$self->getMS132()."','".$self->getMS3()."','".$self->getMnS3()."','".$self->getMnS121()."','".$self->getMnS132()."','".$self->getNxTSqCE()."', '".$self->getMnS152()."');\n";
	return $sql;
}
#".$self->get().",

1;
