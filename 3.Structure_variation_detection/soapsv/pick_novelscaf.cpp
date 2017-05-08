#include<iostream>
#include<fstream>
#include<string>
#include<map>

using namespace std;

map<string, int> data;
map<string, int>::iterator data_iter;

const char* README=
"argv1: all_scaflist\n"
"argv2: selected scaflist\n";

int main(int argc, char** argv)
{
	if(argc<3)
	{
		cerr<<README;
		exit(1);
	}
	
	for(;;)
	{
		string stemp;
		ifstream I(argv[1]);
		for(;;)
		{
			getline(I, stemp, '\n');
			if(!I) break;
		//	stemp = stemp.substr(0,string::npos); // stemp.find(" ")-1);
			data.insert(make_pair<string, int>(stemp, 0));
		}
		break;
	}

	for(;;)
	{
		string stemp;
		ifstream I(argv[2]);
		for(;;)
		{
			getline(I, stemp, '\n');
			if(!I) break;
			if(data.find(stemp) == data.end())
				cerr<<stemp<<" not found!"<<endl;
			++data[stemp];
		}
		break;
	}

	for(data_iter = data.begin(); data_iter != data.end(); ++ data_iter)
	{
		if(data_iter->second == 0)
			cout<<data_iter->first<<endl;
	}
}
