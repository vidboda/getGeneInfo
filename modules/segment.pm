package segment;
sub strict;

sub new {
	my ($class, $nm, $gene) = @_;
        my $self;
	$self = {
                'NM' => $nm,
		'geneName' => $gene,
		'type' => 'intron',
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
sub setType {my ($self, $type) = @_;$self->{type} = $type}
sub getType {my $self = shift;return $self->{type}}
sub setNumber {my ($self, $num) = @_;$self->{num} = $num;$self->{name} = $num;}#default name = number
sub getNumber {my $self = shift;return $self->{num}}
sub setName {my ($self, $name) = @_;$self->{name} = $name}
sub getName {my $self = shift;return $self->{name}}
sub setSize {my ($self, $size) = @_;$self->{size} = $size}
sub getSize {my $self = shift;return $self->{size}}
sub setStartPhase {my ($self, $start_phase) = @_;$self->{start_phase} = $start_phase}
sub getStartPhase {my $self = shift;return $self->{start_phase}}
sub setEndPhase {my ($self, $end_phase) = @_;$self->{end_phase} = $end_phase}
sub getEndPhase {my $self = shift;return $self->{end_phase}}
sub setExonFrame {my ($self, $exon_frame) = @_;if ($exon_frame == '0'){$self->{exon_frame} = 'zero'}else {$self->{exon_frame} = $exon_frame}}
sub getExonFrame {my $self = shift;return $self->{exon_frame}}
sub setStartG {my ($self, $start_g) = @_;$self->{start_g} = $start_g}
sub getStartG {my $self = shift;return $self->{start_g}}
sub setEndG {my ($self, $end_g) = @_;$self->{end_g} = $end_g}
sub getEndG {my $self = shift;return $self->{end_g}}
sub setStartG38 {my ($self, $start_g38) = @_;$self->{start_g38} = $start_g38}
sub getStartG38 {my $self = shift;return $self->{start_g38}}
sub setEndG38 {my ($self, $end_g38) = @_;$self->{end_g38} = $end_g38}
sub getEndG38 {my $self = shift;return $self->{end_g38}}
#sub set {my ($self, ) = @_;$self->{} = }
#sub get {my $self = shift;return $self->{}}

sub toPrint {
	my $self = shift;
	my $exon_frame;
	if ( $self->getExonFrame()) {
		$exon_frame = $self->getExonFrame();
		if ($exon_frame eq 'zero') {$exon_frame = 0}
	}
	my $txt = $self->getGeneName()."\t".$self->getNM()."\t".$self->getType()."\t".$self->getNumber()."\t".$self->getName()."\t".$self->getSize()."\t".$self->getStartG()."\t".$self->getEndG();
	if ($self->getStartG38()) {$txt .= "\t".$self->getStartG38()."\t".$self->getEndG38()}
	if (defined($exon_frame)) {$txt .= "\t$exon_frame"}
	$txt .= "\n";
	return $txt;
}


sub toBed {
	my ($self, $chr, $strand, $offset) = @_;
	if ($strand eq '+') {return "$chr\t".($self->getStartG()-$offset)."\t".($self->getEndG()+$offset)."\t".$self->getNM().'-'.$self->getNumber()."\t0\t$strand\n"}
	else {return "$chr\t".($self->getEndG()-$offset)."\t".($self->getStartG()+$offset)."\t".$self->getNM().'-'.$self->getNumber()."\t0\t$strand\n"}
}

sub toSQL {
	my $self = shift;
	my $exon_frame = $self->getExonFrame();
	if ($exon_frame eq 'zero') {$exon_frame = '0'}
	my $sql = "INSERT INTO segment (nom_gene, type, numero, taille, nom, start_g, end_g, start_g_38, end_g_38, exon_frame) VALUES ('{\"".$self->getGeneName()."\",\"".$self->getNM()."\"}','".$self->getType()."','".$self->getNumber()."','".$self->getSize()."','".$self->getName()."','".$self->getStartG()."','".$self->getEndG()."','".$self->getStartG38()."','".$self->getEndG38()."','$exon_frame');\n";
	return $sql;
}
#replace start_phase, end_phase from ensembl with exon_frames from ucsc in u2
#ALTER ALTER segment ADD exon_frame VARCHAR(2);
#ALTER ALTER segment ADD CONSTRAINT exon_frame_check CHECK ((exon_frame)::text = ANY(ARRAY[('-1'::character varying)::text,('0'::character varying)::text,('1'::character varying)::text,('2'::character varying)::text]));


1;