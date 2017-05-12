
# This file processes the configurations options in Config.sh, producing 
# two files:
#
#   gadgetconfig.h         to be included in each source file (via allvars.h)
#   compile_time_info.c    code to be compiled in, which will print the configuration 
#
if(@ARGV[0])
{
   open(FILE, @ARGV[0]);
}
else
{
   open(FILE, "Config.sh");
}
open(OUTFILE, ">gadgetconfig.h");
open(COUTF,   ">compile_time_info.c");

print COUTF "#include <stdio.h>\n";
print COUTF "void output_compile_time_options(void)\n\{\n";
print COUTF "printf(\n";

while($line=<FILE>)
{
    chop $line;

    @fields = split ' ' , $line;

    if(substr($fields[0], 0, 1) ne "#")
    {
	if(length($fields[0]) > 0)
	{
	    @subfields = split '=', $fields[0];

	    print OUTFILE "#define $subfields[0] $subfields[1]\n";
            print COUTF   "\"        $fields[0]\\n\"\n";
	}
    }
}

print COUTF "\"\\n\");\n";
print COUTF "\}\n";
