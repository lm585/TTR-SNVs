/*

100     Chr1    21089763        21089921        Chr1    21090056        21090218        3       F,      R,      INSERTION       380.333333333333        UNBAL   478    13       10      101,204,3,34,41,48,89,96,97,G13,

output 
Chr1 21089921 21090056  INSERTION     0	1	0 0 0 0 0 1 1 1 0 1 0 ....

126 accessions, if presented 1; otherwise 0;


vector<string> accessionVect;      //input file containing 126 accession names
construct map<string, bool> acc;
output header
for each line in conserv.txt file  //representative conserv.txt file
  parse line into fields
  call getAccPresent()
  output selected fields
  for each elem in accessionVector //output file
    if(acc[elem] == true)
      print "\t" 1
    else
      print "\t" 0


getAccPresent( map<string, bool> & acc, string col1, string field)
{
 set each member of acc to be false
 if(acc.count(col1) == 0)  //a new accession not in the input file
   cerr << new accession 
   exit(1)
 acc[col1] = 'true';
 parse field into a list
 for each member m in the list
   if(acc.count(col1) == 0)
     cerr << new accession
     exit(1)
   acc[m] = 'true'
}
 */

#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <set>
#include <vector>
#include <sstream>
#include <map>

using namespace std;

struct chrSeq
{
 string name;
 vector<char> seq;
};

vector<chrSeq> glob_chr;
int glob_lineNum;


void setChrSeq(ifstream & fasta);
void readFile(const char * argv, ifstream & input);
void writeFile(const char * argv, ofstream & output);
void getFieldContent(vector<string> & fields, char del, string line);
void getFieldContent2(vector<string> & fields, string del, string line);

int main(int argc, char *argv[])
{
 if(argc != 4)
 {
  cerr << "Usage: executable promoterSeqFile geneNameFile outputFile" << endl;
  cerr << "all seq will be capitolized. one sequence one line" << endl;
  return 1;
 }

 ifstream input;
 ofstream output;
 string line;
 vector<string> lineFields;
 map<string, string> geneID;

 readFile(argv[1], input);
 setChrSeq(input);
 input.close();
 //for(int i = 0; i < glob_chr.size(); i++)
 {
  int i = glob_chr.size() - 1;
  cout << glob_chr.size() << '\t' << glob_chr[i].name << '\t' << glob_chr[i].seq.size() << endl;
 }
 
 readFile(argv[2], input);
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent2(lineFields, " \t\n", line);
  for(int i = 0; i < lineFields.size(); i++)
  {
   geneID[lineFields[i]] = "0";
  }
  getline(input, line);
 }
 input.close();
 
 writeFile(argv[3], output);
 for(int i = 0; i < glob_chr.size(); i++)
 {
  string str;
  str = glob_chr[i].name;
  if(geneID.count(str) > 0) //the gene/chr name is in the user-wanted file
  {
   geneID[str] = "1";
   output << '>' << str << endl;
   for(int j = 0; j < glob_chr[i].seq.size(); j++)
   {
    output << (char) toupper(glob_chr[i].seq[j]); 
   }
   output << endl;
  }
 }
 output.close();
 for(map<string, string>::iterator it = geneID.begin(); it != geneID.end(); it++)
 {
  if(it -> second == "0")
  {
   cout << it -> first << " not in the fasta file" << endl;
  }
 }
 return 0;
}

void setChrSeq(ifstream & fasta)
{
 string line, str;
 char * unitStr, * lineCstr;
 vector<string> lineStrList;
 chrSeq myChrSeq;

 getline(fasta, line);
 while(!fasta.eof())
 {
  lineStrList.clear();
  lineCstr = new char [line.length() + 10];
  strcpy(lineCstr, line.c_str());
  unitStr = strtok(lineCstr, " \t\n");
  while(unitStr != NULL)
  {
   str = unitStr;
   lineStrList.push_back(str);
   unitStr = strtok(NULL, " \t\n");
  }

  if(lineStrList.size() > 0)
  {
   if(lineStrList[0][0] == '>')
   {
    if(myChrSeq.seq.size() > 0)
    {
      glob_chr.push_back(myChrSeq);
    }
    myChrSeq.name = lineStrList[0].substr(1);
    myChrSeq.seq.clear();
   }
   else
   {
    for(int i = 0; i < lineStrList.size(); i++)
      for(int j = 0; j < lineStrList[i].length(); j++)
        myChrSeq.seq.push_back(lineStrList[i][j]);
   }
  }
  getline(fasta, line);
 }
 if(myChrSeq.seq.size() > 0)
 {
  glob_chr.push_back(myChrSeq);
 }
 return;
}

void readFile(const char * argv, ifstream & input)
{
 input.open(argv, ios::in);
 if(!input)
 {
  cerr << argv << " can not be opened for reading.\n";
  exit(1);
 }
}

void writeFile(const char * argv, ofstream & output)
{
 output.open(argv, ios::out);
 if(!output)
 {
  cerr << argv << " can not be opened for writing.\n";
  exit(1);
 }
}

void getFieldContent(vector<string> & fields, char del, string line)
{
 vector<int> pos;
 string str;
 int len, size;

 fields.clear();
 for(int i = 0; i < line.length(); i++)
 {
  if(del == line[i])
    pos.push_back(i);
 }
 if(pos.size() > 0)
 {
  len = pos[0];
  str = line.substr(0, len);
  if(len > 0)
    fields.push_back(str);
  for(int i = 0; i < pos.size() - 1; i++)
  {
   len = pos[i+1] - pos[i] - 1;
   str = line.substr(pos[i] + 1, len);
   fields.push_back(str);
  }
  size = pos.size();
  if(pos[size-1] < line.length() - 1) //not at the end of line
  {
   str = line.substr(pos[size-1] + 1);
   fields.push_back(str);
  }
 }
 else
 {
  fields.push_back(line);
 }
}

void getFieldContent2(vector<string> & fields, string del, string line)
{
 vector<int> pos;
 string str;
 int len, size;

 fields.clear();
 for(int i = 0; i < line.length(); i++)
 {
  if(del.find(line[i]) != string::npos)
    pos.push_back(i);
 }
 if(pos.size() > 0)
 {
  len = pos[0];
  str = line.substr(0, len);
  if(len > 0)
    fields.push_back(str);
  for(int i = 0; i < pos.size() - 1; i++)
  {
   len = pos[i+1] - pos[i] - 1;
   if(len > 0)
   {
    str = line.substr(pos[i] + 1, len);
    fields.push_back(str);
   }
  }
  size = pos.size();
  if(pos[size-1] < line.length() - 1) //not at the end of line
  {
   str = line.substr(pos[size-1] + 1);
   fields.push_back(str);
  }
 }
 else
 {
  fields.push_back(line);
 }
}


