/**
    bd-anchors: Bidirectional String Anchors
    Copyright (C) 2021 Grigorios Loukides and Solon P. Pissis. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <iostream>
#include <unordered_map>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include <unordered_set>
#include "utils.h"
using namespace std;

/* Booth's O(n)-time algorithm -- slightly adapted for efficiency */
INT minlexrot( string &X, INT *f, INT n )
{  
	INT n_d = n<<1;
  	for(INT i = 0; i < n_d; ++i)	f[i] = (INT) -1;

  	INT k = 0;
  	for (INT j = 1; j < n_d; ++j)
  	{
                unsigned char sj = X[j%n];
                INT i = f[j - k - 1];
                while (i != (INT)-1 && sj != X[(k + i + 1)%n])
                {
                        if (sj < X[(k + i + 1)%n])        k = j - i - 1;
                        i = f[i];
                }
				
                if (i == (INT) - 1 && sj != X[(k + i + 1)%n])
                {
                        if (sj < X[(k+i+1)%n])    k = j;
                        f[j - k] = -1;
                }
                else
                        f[j - k] = i + 1;
   	}
   	return k;
}

/* Computes the bd-anchors of a string of length n in O(n.ell) time */
void fast_anchors(string &whole_string, unordered_set<INT> &my_map, INT ell)
{
	INT whole_string_len = whole_string.size();
   	INT start_pos = 0;

        INT *f = new INT[ell<<1];
	while( (start_pos + ell) <= whole_string_len )
  	{
		string string_in_win = whole_string.substr(start_pos, ell);
		INT anchor_pos = start_pos + minlexrot(string_in_win, f, ell);
	 	my_map.insert(anchor_pos);
     		start_pos = start_pos + 1;
  	}
	delete []f;
}

/* Booth's O(n)-time algorithm -- slightly adapted for efficiency -- neglecting the last r rotations */
INT red_minlexrot( string &X, INT *f, INT n, INT r )
{  
	INT n_d = n<<1;
  	for(INT i = 0; i < n_d; ++i)	f[i] = (INT) -1;

  	INT k = 0;
  	for (INT j = 1; j < n_d; ++j)
  	{
                unsigned char sj = X[j%n];
                INT i = f[j - k - 1];
                while (i != (INT)-1 && sj != X[(k + i + 1)%n])
                {
                        if (sj < X[(k + i + 1)%n] && j - i - 1 < n - r )        k = j - i - 1;
                        i = f[i];
                }
				
                if (i == (INT) - 1 && sj != X[(k + i + 1)%n])
                {
                        if (sj < X[(k+i+1)%n] && j - i - 1 < n - r )    k = j;
                        f[j - k] = -1;
                }
                else
                        f[j - k] = i + 1;
   	}
   	return k;
}

/* Computes the reduced bd-anchors of a string of length n in O(n.ell) time */
void red_fast_anchors(string &whole_string, unordered_set<INT> &my_map, INT ell, INT r)
{
	INT whole_string_len = whole_string.size();
   	INT start_pos = 0;

        INT *f = new INT[(ell<<1)];
	while( (start_pos + ell) <= whole_string_len )
  	{
		string string_in_win = whole_string.substr(start_pos, ell);
		INT anchor_pos = start_pos + red_minlexrot(string_in_win, f, ell, r);
	 	my_map.insert(anchor_pos);
     		start_pos = start_pos + 1;
  	}
	delete []f;
}

/* Computes the size of the input alphabet in linear time */
INT alphabet_size(string s)
{
    unordered_map<unsigned char, INT> h;
    for (INT i = 0; i < s.length(); i++)	h[s[i]]++;
    return h.size();
}

int main(int argc, char **argv)
{
	if( argc != 3 )
 	{
        	cout<<"Wrong arguments!\n";
 		cout<<"./construct <text_file> <ell>\n";
 		exit(-1);
 	}

 	ifstream is;
 	is.open (argv[1], ios::in | ios::binary);

 	std::string str2(argv[2]);
 	INT ell;
 	std::stringstream(str2)>>ell;
 	

 	cout<<"The parameter ell is set to "<<ell<<endl;
	
	ifstream in_file(argv[1], ios::binary);
   	in_file.seekg(0, ios::end);
   	INT file_size = in_file.tellg();

  	vector<unsigned char> text;
  	char c = 0;
	for (INT i = 1; i < file_size; i++)
	{
		is.read(reinterpret_cast<char*>(&c), 1);
		text.push_back( (unsigned char) c );
	}
  	is.close();

	string text_string(text.begin(), text.end());
	INT sigma = alphabet_size(text_string);
 	cout<<"The input alphabet is of size "<<sigma<<endl;

	if( ell < 1 || ell > text_string.size())
      	{
        	cout<<"ell must be in [1, "<<text_string.size()<<"]"<<endl;
         	return -1;
      	}

 	unordered_set<INT> text_anchors;
      	fast_anchors(text_string, text_anchors, ell);
	INT g1 = text_anchors.size();
	cout<<"The text is of length "<< text_string.size() << " and has "<<g1<<" bd-anchors of order "<<ell<<endl;
	cout<<"The density is "<<(double) g1 / text_string.size()<<endl;
 	
 	INT r = ceil(3*log2(ell)/log2(sigma));
	if( r < 0 || r > ell - 1)
      	{
        	cout<<"The parameter r must be in [0,ell-1]"<<endl;
         	r = 0;
      	}
 	cout<<"The parameter r is set to "<<r<<endl;
	unordered_set<INT> text_red_anchors;
	red_fast_anchors(text_string, text_red_anchors, ell, r);
	INT g2 = text_red_anchors.size();
	cout<<"The text is of length "<< text_string.size() << " and has "<<g2<<" reduced bd-anchors of order "<<ell<<endl;
	cout<<"The reduced density is "<<(double) g2 / text_string.size()<<endl;

	return 0;
}
