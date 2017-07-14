/*
input file 1:
NP_000362.1     Human
NP_001009137.1  Chimp
XP_001099005.1  monkey
XP_537290.1     dog
NP_776392.1     cattle


input file 2:
CLUSTAL multiple sequence alignment by MUSCLE (3.8)


NM_013697.5         CTAATCTCCCTAGGCAAGGTTCATATTTGTGTAGGTTACTTATTCTCCTTTTGTTGACTA
NM_012681.2         ------------------------------------------------------------
XM_002713486.3      ------------------------------------------------------------
NM_001261679.1      ------------------------------------------------------------
NM_000371.3         ----------------------------------------------------GTTGACTA
NM_001009137.1      ------------------------------------------------------------
XM_001495182.3      ------------------------------------------------------------
NM_173967.3         ------------------------------------------------------------
NM_001009800.1      ------------------------------------------------------------
NM_214212.1         ------------------------------------------------------------
XM_537290.4         ------------------------------------------------------------
NM_001287556.1      ------------------------------------------------------------


NM_013697.5         CCTGACAGGATGGCTTCCCTTCGACTCTTCCTCCTTTGCCTCGCTGGACTGGTATTTGTG
NM_012681.2         CCTGACAGGATGGCTTCCCTTCGCCTGTTCCTCCTCTGCCTCGCTGGACTGATATTTGCG
XM_002713486.3      CTTGA-AGCATGGCTTCTCAGCACCTCCTCCTCCTTTGCCTCGCTGGACTGGTATTTGTG
NM_001261679.1      CTTGGCAGGATGGCTTCTCATCGTCTGCTCCTCCTCTGCCTCGCTGGACTGGTATTTGTG
NM_000371.3         CTTGGCAGGATGGCTTCTCATCGTCTGCTCCTCCTCTGCCTTGCTGGACTGGTATTTGTG
NM_001009137.1      CTTGGCAGGATGGCTTCTCATCGTCTGCTTCTCCTCTGCCTTGCTGGACTGGTATTTGTG
XM_001495182.3      CTTGGCACGATGGCTTCTCATCACCTGCTTCTCCTCTGCCTCGCCGGACTGGTATTTGCA
NM_173967.3         CTTGGCAGGATGGCTTCCTTCCGTCTGTTCCTCCTTTGTCTCGCTGGACTGGTGTTTGTG
NM_001009800.1      CTTGGCAGGATGGCTTCCTTCCGTCTGCTCCTCCTTTGCCTCGCTGGACTGGTGTTTGTG
NM_214212.1         CTTGGTAGGATGGCTTCTTACCGTCTGCTCCTCCTTTGCCTCGCTGGACTGGTGTTTGTG
XM_537290.4         CTTGGCAGGATGGCTTCTCATCGTCTGCTCCTCCTCTGCCTCGCTGGACTGGTACTTGTG
NM_001287556.1      CTTGGCAGGATGGCTTCTCA---TCTGCTCCTCCTTTGCCTTGCTGGACTGGTACTTGTG
                    * **  *  ********       **  * ***** ** ** ** ****** *  ***

NM_013697.5         TCTGAAGCTGGCCCCGCGGGTGCTGGAGAATCCAAATGTCCTCTGATGGTCAAAGTCCTG
....

input:
10(length of "human     ")
42(amino acid position)
5(extend left 5 positions, extend right 5 positions)
output-file-name


output:
taxa output order same as in input file 1
human     LDAVRGSPAIN-X
Chimp     ...........-.
monkey    .........V.-.
dog       .........V.-.
cattle    .........A.H-


struct row{
  vector<char> cell;
};
vector<row> matrix;
vector<string> matrixTaxons;
vector<int> aaPos; 

//matrix[0] MASHRLLLLCLAGLVFVSEAGPT----GT....E--
//matrix[0].cell[0] M; cell[1] A; cell[2] S; ...., cell[150] E; cell[151] -; cell[152] -;
//aaPOs[0] 1 //human AA start at pos 1
//aaPOs[1] 2
//aaPOs[150] 147; //4 gaps inserted in the human aa 
//aaPOs[151] 147;
//aaPOs[152] 147;
//to find the alignment block for human aa 20..25..30
//scan aaPos[] : aaPos[19] == 20; startIndex = 19;
//aaPos[33]==30; endIndex = 33;
//extract matrix[0].cell[19] ..matrix[0].cell[33] //15-bp
//        matrix[1].cell[19] ..matrix[1].cell[33]
//        matrix[y].cell[19] ..matrix[y].cell[33]
//        for output and conservation statistics
map<string, string> mapProteinTaxon; 

input-file1 : human - chimp - monkey - dog  (vectTaxonOrder)
multiple align block order:  monkey - dog - human -chimp  (matrixTaxons)
                             matrix[0] - [1] - [2] - [3]

output: matrix[2] - [3] - [0] - [1]
   matrixTaxons[2] - [3] - [0] - [1]
design vector<int> opTaxonIndex;
for taxon in vectTaxonOrder
  find the index in matrixTaxons
  if found, push_back the index to opTaxonIndex

human 2
chimp 3
monkey 0
dog 1



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

struct row{
  vector<char> cell;
};

struct param{
 int startPos, targetPos, endPos; 
 string anchorProteinID;
};

void readFile(const char * argv, ifstream & input);
void writeFile(const char * argv, ofstream & output);
void getFieldContent(vector<string> & fields, char del, string line);
void getFieldContent2(vector<string> & fields, string del, string line);

int main(int argc, char *argv[])
{
 if(argc != 7)
 {
  cerr << argv[0] << " input_file1(see C++ header)  input_file2      output_file    5(the number of spaces between longest taxon-name and alignment [4])		42(amino acid position [5])";
  cerr << "	5(extend left 5 positions, extend right 5 positions [6]) "; 
  cerr <<  endl;
  cerr << "*************************************************************" << endl;
  cerr << " In each multi-alignment block, the order of taxa (Human, Chimp, monkey...) must be same" << endl;
  cerr << "output: replace identical aa with '.' " << endl; 
  return 1;
 }

 ifstream input;
 ofstream output;
 string line;
 vector<string> lineFields;
 map<string, string> mapProteinTaxon;
 vector<row> matrix;
 vector<string> matrixTaxons, vectTaxonOrder;
 vector<int> aaPos, vectConserv, opTaxonIndex;
 set<string> setTaxonF;
 param myParam;
 int k, numOfProt, startIndex, endIndex, numOfSpace ;

 myParam.targetPos = atoi(argv[5]);
 myParam.startPos = myParam.targetPos - atoi(argv[6]);
 myParam.endPos = myParam.targetPos + atoi(argv[6]);
 myParam.anchorProteinID = "";
 readFile(argv[1], input);
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent2(lineFields, "\t\n ", line);
  if(lineFields.size() == 2) 
  {
   mapProteinTaxon[lineFields[0]] = lineFields[1];//ID -> taxon
   vectTaxonOrder.push_back(lineFields[1]);
   if(myParam.anchorProteinID == "")
     myParam.anchorProteinID = lineFields[0]; //assigned to "NM_012345.2"; next iteration, != ""
  }
  getline(input, line);
 }
 input.close();
 
 readFile(argv[2], input);
 getline(input, line);
 while(!input.eof())
 {
  //set<string> setTaxonF;
  string str1, str2;
  getFieldContent2(lineFields, "\t\n ", line);
  if(lineFields.size() == 2 && mapProteinTaxon.count(lineFields[0]) == 1 )
  //two fields; 1st field is an protein ID stored in input file 1
  {
   str1 = lineFields[0];//protein ID of a taxon
   if(mapProteinTaxon.count(str1) == 0) //protein ID has no corresponding taxon
   {
    cerr << str1 << " has no corresponding taxon in the file " << argv[1] << endl;
    return 1;
   }
   matrixTaxons.push_back(mapProteinTaxon[str1]); //store taxon info in the align block order
   cout << "taxon in the alignmnet block" << '\t' << str1 << '\t' << mapProteinTaxon[str1] << endl; 
   if(setTaxonF.find(str1) != setTaxonF.end()) //find a taxon second time; begin of 2nd align block
   {
    numOfProt = setTaxonF.size();
    break;
   }
   else
   {
    setTaxonF.insert(str1);
   }
  }
  getline(input, line);
 }
 input.close();
 // build opTaxonIndex vector
 {
  for(int i = 0; i < vectTaxonOrder.size(); i++)
  {
   for(int j = 0; j < matrixTaxons.size(); j++)
   {
    if(matrixTaxons[j] == vectTaxonOrder[i] )
    {
     opTaxonIndex.push_back(j);
     break;
    }
   }
  }
  cout << "opTaxonIndex size: " << opTaxonIndex.size() << endl;
 }

 row tempRow;
 for(int i = 0; i < numOfProt; i++)
 {
  matrix.push_back(tempRow);
 }
 k = -1;
 readFile(argv[2], input);
 getline(input, line);
 while(!input.eof())
 {
  string str1, str2;
  int pos;
  getFieldContent2(lineFields, "\t\n ", line);
  if(lineFields.size() == 2 && mapProteinTaxon.count(lineFields[0]) == 1)
  {
   str1 = lineFields[1]; //VSEAGPT----GTGESK
   //pos = 1; //atoi(lineFields[1].c_str());
   k++;
   for(int i = 0; i < str1.length(); i++)
   {
    matrix[k % numOfProt].cell.push_back(str1[i]);
    if(lineFields[0] == myParam.anchorProteinID)
    {
     if(aaPos.size() == 0)
     {
      if(str1[0] == '-')
        aaPos.push_back(0);
      else
        aaPos.push_back(1);
     }
     else
     {
      int j = aaPos.size() - 1; //size 6 , [5]
      int t = aaPos[j];
      if(str1[i] == '-')
        aaPos.push_back(t);
      else
        aaPos.push_back(t+1);
     }
    }
   }
  }
  getline(input, line);
 }
 input.close();
 for(int i = 0; i < numOfProt; i++)
 {
  cout << matrixTaxons[i] << '\t' << matrix[i].cell.size() << endl;
 }
 cout << "aaPos.size() =" << aaPos.size() << endl;

 if(myParam.startPos < aaPos[0]) myParam.startPos = aaPos[0];
 k = aaPos.size() - 1;
 if(myParam.endPos > aaPos[k]) myParam.endPos = aaPos[k];
 startIndex = -1;
 endIndex = -1;
 for(int i = 0; i < aaPos.size(); i++)
 {
  if(aaPos[i] == myParam.startPos && startIndex == -1)
    startIndex = i;
  if(aaPos[i] == myParam.endPos && endIndex == -1)
    endIndex = i;
 }
//vectConserv
 k = opTaxonIndex[0];
 for(int i = startIndex; i <= endIndex; i++)
 {
  vectConserv.push_back(0);
 }
 for(int i = 0; i < numOfProt; i++)
 {
  for(int j = startIndex; j <= endIndex; j++)
  {
   if(toupper(matrix[i].cell[j]) == toupper(matrix[k].cell[j]))
   {
    vectConserv[j-startIndex]++;
   }
  }
 }
 writeFile(argv[3], output);
 //numOfSpace
 k = 0;
 for(int i = 0; i < matrixTaxons.size(); i++)
 {
  if(matrixTaxons[i].length() > k)   k = matrixTaxons[i].length();
 }
 numOfSpace = k + atoi(argv[4]);
 for(int i = 0; i < numOfSpace; i++)
 {
  output << ' ';
 }
 for(int j = startIndex; j <= endIndex; j++)
 {
  output << vectConserv[j-startIndex] << " ";
 }
 output << endl;
 for(int i = 0; i < numOfSpace; i++)
 {
  output << ' ';
 }
 for(int j = startIndex; j <= endIndex; j++)
 {
  if(aaPos[j] == myParam.targetPos)
    output << '|';
  else
    output << ' ';
 }
 output << endl;
 for(int i = 0; i < opTaxonIndex.size(); i++)
 {
  k = opTaxonIndex[i]; //k in: 2 - 3 - 0 - 1
  int anch = opTaxonIndex[0]; //anch always is 2
  output << matrixTaxons[k];
  for(int t = 0; t < numOfSpace - matrixTaxons[k].length(); t++)
  {
   output << ' ';
  }
  for(int j = startIndex; j <= endIndex; j++)
  {
   if(i == 0)
     output << matrix[anch].cell[j]; 
   else
   {
    if(matrix[anch].cell[j] == '-')
      output << matrix[k].cell[j];
    else
    {
     if(toupper(matrix[k].cell[j]) == toupper(matrix[anch].cell[j]))
       output << '.';
     else
       output << matrix[k].cell[j];
    }
   }
  }
  output << endl;
 }

 output.close();
 return 0;
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


