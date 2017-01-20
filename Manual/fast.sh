###Non-looping version
#Insert in the line before "close SHL"
	print SHL "echo You may proceed! >> FLAGFILE.txt";

#Insert after "system qsub ... "
	while (not -e "FLAGFILE.txt") {sleep 5};
	system("rm -f FLAGFILE.txt");
	
	
#### Loop version
#Before loop
system("rm -f flag*");

#Last line before newline SHL print or close
print SHL "echo You may proceed! >> flag$sub.txt";
#Check for if-else statements at loop end

#End of script
foreach my $sub (@subs) {
while (not -e "flag$sub.txt") {sleep 5};
}