package domain;
use strict;

sub new {
	my ($class, $nm, $gene) = @_;
        my $self;
	$self = {
            "NM"=> $nm,
			'geneName' => $gene
		};
	bless ($self, $class);
	return $self;
}

sub setNM {
        my ($self, $nm) = @_;
        $self->{NM} = $nm;
}
sub getNM {
        my $self = shift;
        return $self->{NM};
}
sub setGeneName{my ($self, $gene) = @_;$self->{geneName} = $gene}
sub getGeneName{my $self = shift;return $self->{geneName}}
sub setName {my ($self, $name) = @_;$self->{name} = $name}
sub getName {my $self = shift;return $self->{name}}
sub setStartAA {my ($self, $start_aa) = @_;$self->{start_aa} = $start_aa}
sub getStartAA {my $self = shift;return $self->{start_aa}}
sub setEndAA {my ($self, $end_aa) = @_;$self->{end_aa} = $end_aa}
sub getEndAA {my $self = shift;return $self->{end_aa}}
#sub set {my ($self, ) = @_;$self->{} = }
#sub get {my $self = shift;return $self->{}}

sub toPrint {
	my $self = shift;
	return $self->getGeneName()."\t".$self->getNM()."\t".$self->getName()."\t".$self->getStartAA()."\t".$self->getEndAA()."\n";
}

sub toLOVD {
	my $self = shift;
	return $self->getName()." (".$self->getStartAA()."-".$self->getEndAA().")\n";
}

sub toSQL {
	my ($self, $short_prot, $system) = @_;
	my $sql = '';
	if ($system eq 'md') {
		$sql = "INSERT INTO protein_domain (name, gene_name, aa_start, aa_end) VALUES ('".$self->getName()."','{\"".$self->getGeneName()."\",\"".$self->getNM()."\"}','".$self->getStartAA()."','".$self->getEndAA()."');\n"
	}
	elsif ($system eq 'u2') {
		$sql = "INSERT INTO domaine (nom, nom_prot, aa_deb, aa_fin) VALUES ('".$self->getName()."','".$short_prot."','".$self->getStartAA()."','".$self->getEndAA()."');\n"
	}
	return $sql;
}

1;