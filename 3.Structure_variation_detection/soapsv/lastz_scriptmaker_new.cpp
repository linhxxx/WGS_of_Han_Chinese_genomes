#include<iostream>
#include<fstream>
#include<string>
using namespace std;

const char* README=
"argv1 fa list\n"
"argv2 output dir\n";

int main(int argc, char** argv)
{
	if(argc==1)
	{
		cerr<<README;
		exit(1);
	}
	ifstream I_fa(argv[1]);
	string sf;
	for(;;)
	{
		I_fa>>sf;
		if(!I_fa) break;
		string sf_sub,output;
		sf_sub = sf.substr(sf.rfind("/")+1, string::npos);
		output = sf_sub.substr(0,sf_sub.find("."));
		sf_sub = sf_sub.substr(0,sf_sub.find("_"));
		cout<<"/ifs2/POPULATION/GROUP/zhanghk/project/soft/scr/lastz "	
			<<"--targetcapsule="<<"/ifs2/BC_IP/PROJECT/zhk_test/sv/prj_sv_renal/ref/capsule.list"<<sf_sub<<".fa.capsult_12of19_seed5 "
			<<sf<<" "
			<<"--strand=both --chain --ambiguousn --gapped "
			<<"--ydrop=20000 --gap=1000,1 "
			<<"--noentropy "
			<<"--format=axt "
			<<"--output="<<argv[2]<<"/"<<sf_sub<<"/"<<output<<".axt "
			<<"--markend "
			<<endl;		
	}	
}

