#!/usr/bin/perl

while (my $line = <>)
{
	if ($line =~ /^(>\S+)_(\d+)_\d+\s*(.*)$/) { $header = $1; $st = $2; $desc = $3 }
	elsif ($line =~ /^(>\S+)\s*(.*)$/) { $header = $1; $st = 1; $desc = $2 }
	else
	{
		while ($line =~ /([a-z]+)/g)
		{
			print "$header\_" .($st + pos($line) - length($1))."_".($st + pos($line) - 1)." $desc\n$1\n";
		}
	}
}