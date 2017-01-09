#Insert in the line before "close SHL"
	print SHL "echo You may proceed! >> FLAGFILE.txt";

#Insert after "system qsub ... "
	while (not -e "FLAGFILE.txt") {sleep 5};
	system("rm -f FLAGFILE.txt");