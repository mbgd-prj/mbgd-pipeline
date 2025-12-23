package DomRefine::Draw;
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(draw_tick draw_alignment create_png create_image change_coordinate draw_annotation draw_rect_filled create_colors);

use strict;
use GD;
use DomRefine::Read;
use DomRefine::Motif;
use DomRefine::General;

my $HIGH = 0.7;
my $MIDDLE = 0.5;
my $LOW = 0.3;

sub create_image {
    my ($r_image, $r_color, $a_n, $m, $add, $h_param) = @_;

    ${$h_param}{MARGIN} = 20;
    # ${$h_param}{X_MARGIN} = 135;
    ${$h_param}{X_MARGIN} = 180;
    # ${$h_param}{X_MARGIN} = ${$h_param}{MARGIN} + $gene_length * 7;
    ${$h_param}{X_UNIT} = 1;
    ${$h_param}{Y_UNIT} = 13;
    ${$h_param}{Y_DOMAIN} = 1;
    ${$h_param}{Y_AMINO_ACID} = int(${$h_param}{Y_UNIT} / 2);
    ${$h_param}{Y_CONSERVATION} = ${$h_param}{Y_AMINO_ACID} - 2;

    my $width = ${$h_param}{X_MARGIN} + $m + ${$h_param}{X_MARGIN} + 500;
    $width = max($width, 1200);
    my $height = ${$h_param}{MARGIN} * (@{$a_n} + 1.5) + ${$h_param}{Y_UNIT} * (sum(@{$a_n}) + @{$a_n} * 3 + $add);

    ${$r_image} = new GD::Image($width, $height);
    ${$r_image}->interlaced('true');
    create_colors($r_image, $r_color);
}

sub create_png {
    my ($r_image, $name, %opt) = @_;

    open(PNG, ">$name.png") || die;
    binmode PNG;
    print PNG ${$r_image}->png;
    close(PNG);
    if ($opt{draw_now}) {
	system "display $name.png";
    }
}

sub draw_alignment {
    my ($r_gene, $r_a, $r_b, $r_d, $r_p, $r_get_j, $r_image, $h_color, $r_offset, $h_domain, %opt) = @_;
    my $r_scores_for_line = $opt{r_scores_for_line};
    
    my $N = @{$r_a};
    my $M = @{${$r_a}[0]};

    print STDERR " draw alignment..\n";
    # conservation
    if ($opt{conserv}) {
	for (my $i=0; $i<$N; $i++) {
	    for (my $j=0; $j<$M; $j++) {
		my $color;
		if (${$r_b}[$i][$j] == 1) {
		    if (${$r_p}{$j}{${$r_a}[$i][$j]}>= $HIGH) {
			$color = ${$h_color}{RED};
		    } elsif (${$r_p}{$j}{${$r_a}[$i][$j]} >= $MIDDLE) {
			$color = ${$h_color}{YELLOW};
		    } elsif (${$r_p}{$j}{${$r_a}[$i][$j]} >= $LOW) {
			$color = ${$h_color}{CYAN};
		    }
		}
		if ($color) {
		    draw_tick($r_image, $color, $opt{X_UNIT}*$j, ${$r_offset}+$opt{Y_UNIT}*$i, %opt, y_pos_shift => $opt{Y_CONSERVATION});
		}
	    }
	}
    }
    # amino acid
    for (my $i=0; $i<$N; $i++) {
	for (my $j=0; $j<$M; $j++) {
	    if (${$r_b}[$i][$j] == 1) {
		draw_tick($r_image, ${$h_color}{BLACK}, $opt{X_UNIT}*$j, ${$r_offset}+$opt{Y_UNIT}*$i, %opt, y_pos_shift => $opt{Y_AMINO_ACID});
	    }
	}
    }
    # gene labels
    for (my $i=0; $i<$N; $i++) {
	my $gene = ${$r_gene}[$i];
	my @domains = sort {$a<=>$b} keys %{${$h_domain}{$gene}};
	my $domains = join(",", @domains);
	my $label = "$gene($domains)";
	${$r_image}->string(gdSmallFont, $opt{MARGIN}, $opt{MARGIN}+${$r_offset}+$opt{Y_UNIT}*$i, $label, ${$h_color}{BLACK});
    }
    # score for line
    for (my $i=0; $i<$N; $i++) {
	if (defined ${$r_scores_for_line}[$i]) {
	    my $color = ${$h_color}{BLACK};
	    if (${$r_scores_for_line}[$i] <= -0.5) {
		$color = ${$h_color}{RED};
	    } elsif (${$r_scores_for_line}[$i] <= 0) {
		$color = ${$h_color}{BLUE};
	    }
	    ${$r_image}->string(gdSmallFont, $opt{X_MARGIN}+$opt{X_UNIT}*$M+10, $opt{MARGIN}+${$r_offset}+$opt{Y_UNIT}*$i, ${$r_scores_for_line}[$i], $color);
	}
    }
    # gene_descr
    for (my $i=0; $i<$N; $i++) {
	my $gene_descr = get_geneset_descr(${$r_gene}[$i], mysql => $opt{mysql}, r_gene_descr => $opt{r_gene_descr});
	if ($gene_descr) {
	    ${$r_image}->string(gdSmallFont, $opt{X_MARGIN}+$opt{X_UNIT}*$M+10+120, $opt{MARGIN}+${$r_offset}+$opt{Y_UNIT}*$i, $gene_descr, ${$h_color}{BLACK});
	}
    }

    my @hit_motif;
    my ($motif1, $motif2, @other_motif);
    if ($opt{reference_file}) {
	my @hit_start_j;
	my @hit_end_j;
	my @hit_gene_i;
	my @hit_evalue;

	print STDERR " get motif hits..\n";
	get_hit_positions_sparql($r_gene, $r_get_j, $opt{reference_file}, \@hit_motif, \@hit_gene_i, \@hit_start_j, \@hit_end_j, \@hit_evalue, motifs_to_ignore => $opt{motifs_to_ignore});

	# print STDERR " read $opt{reference_file} ..\n";
	# my %reference_cluster = ();
	# my %reference_domain = ();
	# get_dclst_structure($opt{reference_file}, \%reference_cluster, \%reference_domain);
	# for (my $i=0; $i<$N; $i++) {
	#     get_hit_positions(${$r_gene}[$i], \@hit_start_j, \@hit_end_j, \@hit_motif, \@hit_gene_i, \@hit_evalue, ${$r_get_j}[$i], $i, \%reference_domain, %opt);
	# }

	print STDERR " draw matches to references..\n";
	($motif1, $motif2, @other_motif) = sort_hit_motifs(@hit_motif);
	for (my $k=0; $k<@hit_motif; $k++) {
	    my $color;
	    if ($hit_motif[$k] eq $motif1) {
		$color = ${$h_color}{BLUE};
	    } elsif ($hit_motif[$k] eq $motif2) {
		$color = ${$h_color}{RED};
	    } else {
		$color = ${$h_color}{BLACK};
	    }
	    if (defined $hit_end_j[$k]) {
		draw_rect($r_image, $hit_gene_i[$k], $hit_start_j[$k], $hit_end_j[$k], $color, ${$r_offset}, %opt, y_pos_shift => $opt{Y_DOMAIN});
		# print STDERR "$hit_gene_i[$k] $hit_motif[$k] $hit_start_j[$k]\n";
	    }
	}
    }

    ${$r_offset} += $opt{Y_UNIT} * $N + $opt{MARGIN};

    if ($motif1) {
	${$r_offset} += $opt{MARGIN} * 0.25;
	draw_ledgends($r_image, $h_color, $r_offset, $motif1, $motif2, \@hit_motif, $opt{reference_file}, r_other_motif => \@other_motif, %opt);
    }
}

sub draw_rect_filled {
    my ($r_image, $h_color, $x, $y, %opt) = @_;

    my $color;
    if ($opt{color_idx} == 0) {
	$color = ${$h_color}{PALEBLUE};
    } elsif ($opt{color_idx} == 1) {
	$color = ${$h_color}{PINK};
    } else {
	$color = ${$h_color}{LIGHTGRAY};
    }
    ${$r_image}->filledRectangle(change_coordinate(($opt{MARGIN}+$x, $opt{MARGIN}+$y), (0, 1, 42, $opt{Y_UNIT}-1)), $color);
}

sub draw_annotation {
    my ($r_image, $h_color, $x, $y, %opt) = @_;
    
    my $cluster = "";
    if ($opt{cluster}) {
	$cluster = "[$opt{cluster}]    ";
    }

    if ($opt{homcluster}) {
	$cluster = "$opt{homcluster} $cluster";
    }

    my $scores = "";
    if ($opt{scores}) {
	$scores = "$opt{scores}    ";
    }

    my $annotation = $opt{annotation} || "";

    ${$r_image}->string(gdMediumBoldFont, $opt{MARGIN} + $x, $opt{MARGIN} + $y, "$cluster$scores$annotation", ${$h_color}{BLACK});
}

sub draw_ledgends {
    my ($r_image, $h_color, $r_offset, $motif1, $motif2, $r_hit_motif, $file, %opt) = @_;

    my %name = ();
    my %descr = ();
    my @line = `get_motif_annotation @{$r_hit_motif}`;
    chomp(@line);
    for my $line (@line) {
	my @f = split("\t", $line, -1);
	if (@f != 3) {
	    die $line;
	}
	my ($motif, $name, $descr) = @f;
	$name{$motif} = $name;
	$descr{$motif} = $descr;
    }

    my %ref_hits = ();
    get_annotation($file, \%ref_hits);

    my %count_hit_motif = ();
    for my $motif (@{$r_hit_motif}) {
	$count_hit_motif{$motif} ++;
    }

    if ($motif1) {
	${$r_image}->rectangle(change_coordinate(($opt{MARGIN}, ${$r_offset}), (0, 1, 43, $opt{Y_UNIT}-1)), ${$h_color}{BLUE});
	my $descr = $descr{$motif1} || "";
	${$r_image}->string(gdMediumBoldFont, change_coordinate(($opt{MARGIN}, ${$r_offset}), (50, 0)), "$name{$motif1} ($count_hit_motif{$motif1}/$ref_hits{$motif1}) $descr", ${$h_color}{BLACK});
    	${$r_offset} += $opt{Y_UNIT};
    }
    if ($motif2) {
	${$r_image}->rectangle(change_coordinate(($opt{MARGIN}, ${$r_offset}), (0, 1, 43, $opt{Y_UNIT}-1)), ${$h_color}{RED});
	my $descr = $descr{$motif2} || "";
	${$r_image}->string(gdMediumBoldFont, change_coordinate(($opt{MARGIN}, ${$r_offset}), (50, 0)), "$name{$motif2} ($count_hit_motif{$motif2}/$ref_hits{$motif2}) $descr", ${$h_color}{BLACK});
    	${$r_offset} += $opt{Y_UNIT};
    }
    if (@{$opt{r_other_motif}}) {
	${$r_image}->rectangle(change_coordinate(($opt{MARGIN}, ${$r_offset}), (0, 1, 43, $opt{Y_UNIT}-1)), ${$h_color}{BLACK});
	my @other_motif = ();
	my @descr = ();
	for my $other_motif (@{$opt{r_other_motif}}) {
	    push @other_motif, "$name{$other_motif} ($count_hit_motif{$other_motif}/$ref_hits{$other_motif})";
	    if ($descr{$other_motif}) {
	    	push @descr, $descr{$other_motif};
	    }
	}
	my $descr = join(" | ", @descr) || "";
	${$r_image}->string(gdMediumBoldFont, change_coordinate(($opt{MARGIN}, ${$r_offset}), (50, 0)), "@other_motif $descr", ${$h_color}{BLACK});
    	${$r_offset} += $opt{Y_UNIT};
    }
}

sub draw_rect {
    my ($r_image, $i, $hit_start_j, $hit_end_j, $color, $offset, %opt) = @_;

    ${$r_image}->rectangle(change_coordinate($opt{X_MARGIN}, $opt{MARGIN} + $offset + $opt{Y_UNIT} * $i,
					     $hit_start_j, $opt{y_pos_shift},
					     $hit_end_j, $opt{Y_UNIT}-$opt{y_pos_shift}
					     ), $color);
}

sub draw_tick {
    my ($r_image, $color, $x, $y, %opt) = @_;
    
    ${$r_image}->filledRectangle(change_coordinate($opt{X_MARGIN} + $x, $opt{MARGIN} + $y,
						   0, $opt{y_pos_shift},
						   0, $opt{Y_UNIT}-$opt{y_pos_shift}
						   ), $color);
}

sub create_colors {
    my ($r_image, $r_color) = @_;

    ${$r_color}{WHITE} = ${$r_image}->colorAllocate(255,255,255);
    ${$r_color}{BLACK} = ${$r_image}->colorAllocate(0,0,0);       
    ${$r_color}{GRAY} = ${$r_image}->colorAllocate(128,128,128);
    ${$r_color}{LIGHTGRAY} = ${$r_image}->colorAllocate(160,160,160);

    ${$r_color}{RED} = ${$r_image}->colorAllocate(255,0,0);      
    ${$r_color}{BLUE} = ${$r_image}->colorAllocate(0,0,255);

    ${$r_color}{PINK} = ${$r_image}->colorAllocate(255,182,193);
    ${$r_color}{PALEBLUE} = ${$r_image}->colorAllocate(176,196,222);

    ${$r_color}{MAGENTA} = ${$r_image}->colorAllocate(255,0,255);
    ${$r_color}{YELLOW} = ${$r_image}->colorAllocate(255,255,0);
    ${$r_color}{CYAN} = ${$r_image}->colorAllocate(0,255,255);

    ${$r_color}{VIOLET} = ${$r_image}->colorAllocate(238,130,238);
    ${$r_color}{GREEN} = ${$r_image}->colorAllocate(0,255,0);
    ${$r_color}{DEEPGREEN} = ${$r_image}->colorAllocate(34,139,34);
}

sub change_coordinate {
    my ($x_offset, $y_offset, $x1, $y1, $x2, $y2) = @_;

    my @ret = ();

    push @ret, $x_offset + $x1;
    push @ret, $y_offset + $y1;

    if (defined $x2) {
	push @ret, $x_offset + $x2;
    }
    if (defined $y2) {
	push @ret, $y_offset + $y2;
    }
    
    # return ($x1+$x_offset, $y1+$y_offset, $x2+$x_offset, $y2+$y_offset);
    return @ret;
}

1;
