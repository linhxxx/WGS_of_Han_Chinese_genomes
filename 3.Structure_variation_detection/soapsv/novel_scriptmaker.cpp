#include<iostream>
#include<fstream>
#include<string>
using namespace std;

const char* README=
"argv1 capsule list\n"
"argv2 fa list\n"
"argv3 output dir\n";

int main(int argc, char** argv)
{
	if(argc==1)
	{
		cerr<<README;
		exit(1);
	}
	ifstream I_capsule(argv[1]);
	string sc, sf;
	for(;;)
	{
		I_capsule>>sc;
		string sc_sub;
		sc_sub = sc.substr(sc.rfind("/")+1, string::npos);
		sc_sub = sc_sub.substr(0, sc_sub.find("."));
		if(!I_capsule) break;
		ifstream I_fa(argv[2]);
		for(;;)
		{
			I_fa>>sf;
			string sf_sub;
			sf_sub = sf.substr(sf.rfind("/")+1, string::npos);
			sf_sub = sf_sub.substr(0,sf_sub.find("."));
			if(!I_fa) break;
			cout<<""
				<<"/ifs2/POPULATION/GROUP/zhanghk/project/soft/scr/lastz "	
				<<"--targetcapsule="<<sc<<" "
				<<sf<<" "
				<<"--strand=both --chain --ambiguousn --gapped --ydrop=30000 --gap=1000,1 "
				<<"--format=axt "
				<<"--output="<<argv[3]<<"/"<<sc_sub<<"_"<<sf_sub<<".axt "
				<<"--markend "
				<<endl;
		}		
	}	
}

