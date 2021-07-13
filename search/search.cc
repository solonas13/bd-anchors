#define DEBUG false
#include <algorithm>
#include <chrono>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <sys/time.h>
#include <iostream>
#include <cstdint>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <deque>
#include <boost/functional/hash.hpp>
#include <list>
#include <ctime>
#include <random>
#include <algorithm>

#include "int.h"

#ifdef _USE_64
#include <divsufsort64.h>                                         // include header for suffix sort
#endif
#ifdef _USE_32
#include <divsufsort.h>                                         // include header for suffix sort
#endif

#include <sdsl/bit_vectors.hpp>                                   // include header for bit vectors
#include <sdsl/rmq_support.hpp>

using namespace sdsl;
using namespace std;
using namespace boost;

// Both edit distance computation functions TAKEN FROM https://www.geeksforgeeks.org/edit-distance-dp-5/

INT EditDistDP(string str1, string str2)
{
    INT len1 = str1.length();
    INT len2 = str2.length();
 
 	
	INT **DP=new INT*[len1+1];
	DP[0]=new INT[len1+1];
	DP[1]=new INT[len1+1];
	for(INT i=0;i<len1+1;++i)
	{
		DP[0][i]=0;
		DP[1][i]=0;
	}
 
    for (INT i = 0; i <= len1; i++)
        DP[0][i] = i;
 
    for (INT i = 1; i <= len2; i++) {
        for (INT j = 0; j <= len1; j++) {
            if (j == 0)
                DP[i % 2][j] = i;
            else if (str1[j - 1] == str2[i - 1]) {
                DP[i % 2][j] = DP[(i - 1) % 2][j - 1];
            }            
            else {
                DP[i % 2][j] = 1 + min(DP[(i - 1) % 2][j],
                                       min(DP[i % 2][j - 1],
                                           DP[(i - 1) % 2][j - 1]));
            }
        }
    } 
	INT ret=DP[len2 % 2][len1];
	for(INT i = 0; i < 2; ++i) {
		delete [] DP[i];
	}
	delete [] DP;
	return ret;	
}
          
/* Kasai et al algorithm for LCP construction */
INT LCParray ( unsigned char * text, INT n, INT * SA, INT * ISA, INT * LCP )
{
        INT i=0, j=0;

        LCP[0] = 0;
        for ( i = 0; i < n; i++ ) // compute LCP[ISA[i]]
                if ( ISA[i] != 0 )
                {
                        if ( i == 0) j = 0;
                        else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
                        while ( text[i+j] == text[SA[ISA[i]-1]+j] )
                                j++;
                        LCP[ISA[i]] = j;
                }

        return ( 1 );
}

INT minlexrot( string& X ) 
{    
	//string S = X;
    //	S.append(X.begin(), X.end());
    //	INT n = S.size();
    
    //	unsigned char * seq = (unsigned char *) S.c_str();
	
	X.append(X.begin(), X.end());
	INT n=X.size();
	unsigned char * seq = (unsigned char *)X.c_str();
	
    vector <INT> f (n, (INT)-1);
    
    INT k = 0;
    for (INT j = 1; j < n; ++j) 
	{
       	unsigned char sj = seq[j];
       	INT i = f[j - k - 1];
       	while (i != (INT)-1 && sj != seq[k + i + 1]) 
		{
      		if (sj < seq[k + i + 1])	k = j - i - 1;
           		i = f[i];
       	}

       	if (i==(INT)-1 && sj != seq[k + i + 1]) 
		{
       		if (sj < seq[k+i+1])	k = j;
           		f[j - k] = -1;
       	} 
		else 
       		f[j - k] = i + 1;
    }

    	return k;
}

INT minlexrot( string& X, INT *f, INT n)
{  
  INT n_d=n<<1;
  for(INT i=0;i<n_d;++i)
        f[i]=(INT)-1;

  unsigned char * seq = (unsigned char *) X.c_str();
  INT k=0;
  for (INT j = 1; j < n_d; ++j)
  {
                unsigned char sj = seq[j%n];
                INT i = f[j - k - 1];
                while (i != (INT)-1 && sj != seq[(k + i + 1)%n])
                {
                        if (sj < seq[(k + i + 1)%n])        k = j - i - 1;
                        i = f[i];
                }
				
                if (i==(INT)-1 && sj != seq[(k + i + 1)%n])
                {
                        if (sj < seq[(k+i+1)%n])    k = j;
                        f[j - k] = -1;
                }
                else
                        f[j - k] = i + 1;
   }
   return k;

}

INT minlexrotf( unsigned char *X, INT *f, INT n)
{  
  const INT n_d=n<<1;
  for(INT i=0;i<n_d;++i)
        f[i]=(INT)-1;

  INT k=0;
  for (INT j = 1; j < n_d; ++j)
  {
                const unsigned char sj = X[j%n];
				const INT j_minus_1=j-1;
                INT i = f[j_minus_1-k];
                while (i != (INT)-1 && sj != X[(k + i + 1)%n])
                {
                        if (sj < X[(k + i + 1)%n])        k = j_minus_1- i;
                        i = f[i];
                }
				
                if (i==(INT)-1 && sj != X[(k + i + 1)%n])
                {
                        if (sj < X[(k+i+1)%n])    k = j;
                        f[j - k] = -1;
                }
                else
                        f[j - k] = i + 1;
   }
   return k;

}
/* Computes the bd-anchors of a string of length n in O(wn) time */
void fast_anchors(string &whole_string, unordered_set<INT> &my_map, INT w)
{
	INT whole_string_len=whole_string.size();
   	INT start_pos = 0;

    INT *f=new INT[w<<1];
	while( (start_pos + w) <= whole_string_len )
  	{
		string string_in_win = whole_string.substr(start_pos,w);
		INT anchor_pos = start_pos + minlexrot(string_in_win,f, w);
	 	my_map.insert(anchor_pos);
     		start_pos = start_pos + 1;
  	}
	delete []f;
}

	
//Declaration of pattern matching as it is used in fast_anchors_pattern
void pattern_matching ( string & w, string & a, INT * SA, INT * LCP, std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> &rmq, INT n, INT &pair1, INT &pair2 );
void rev_pattern_matching ( string & w, string & a, INT * SA, INT * LCP, std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> &rmq, INT n, INT &pair1, INT &pair2 );

//THE ARGUMENT whole_string is reversed when this function is called
void fast_anchors_pattern_end(string& whole_string, string &pattern, INT ell, INT * LSA, INT * LLCP, std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> &lrmq, INT * RSA, INT * RLCP, 
std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> &rrmq, INT g,  vector<pair<INT,INT> > &pattern_answers, unordered_map<INT,INT> &anchors_and_window_ends_in_pattern,  unordered_map<INT,INT> &invRSA)
{    
  
   INT pattern_len=pattern.size();
   INT whole_string_len=whole_string.size();
   INT start_pos=0;     
   
   INT *f=new INT[ell<<1];
   //start_pos is in pattern
   while((start_pos+ell) <= pattern_len)
  {    
     unsigned char* string_in_win=(unsigned char*)pattern.substr(start_pos,ell).c_str();  
	//anchor_pos is in pattern
	  //leftmost 
     INT anchor_pos= start_pos+minlexrotf(string_in_win,f,ell);
	
	 auto it2=anchors_and_window_ends_in_pattern.find(anchor_pos);
	 if(it2==anchors_and_window_ends_in_pattern.end())
	 {
		 anchors_and_window_ends_in_pattern.insert(make_pair(anchor_pos,start_pos+ell-1));
	 }
	 else
	 {
		 start_pos++;
		 continue;		 
	 }
	 
	if ( ell+start_pos - anchor_pos >= anchor_pos-start_pos ) //if the right part is bigger than the left part then search the right part to get a smaller interval on RSA
		{
	
  			string right_pattern = pattern.substr(anchor_pos,start_pos+ell-1-anchor_pos+1);
			INT pair1=0,pair2=0;
			pattern_matching ( right_pattern, whole_string, RSA, RLCP, rrmq, g, pair1, pair2 );
  			
			if(pair1>pair2)
			{
				start_pos = start_pos + 1; 
				continue;
			}
			INT anchor_pos_to_return=anchor_pos;
			for(INT i = pair1; i <= pair2; i++ ) //this can be an interval of size 10^10 and only one occurrence is valid.
			{
				
				INT index = RSA[i];				
				anchor_pos=anchor_pos_to_return;
				index--;	anchor_pos--;
				while ( ( anchor_pos >= start_pos ) && ( index >= 0 ) && ( whole_string[index] == pattern[anchor_pos] ) )
				{
					index--; anchor_pos--;
				}
				if(anchor_pos==start_pos-1)
				{
					
					pair<INT,INT> my_pair(anchor_pos_to_return,i);
					pattern_answers.push_back(my_pair);
				}
				
			}
		}
		else //otherise search the left part to get a smaller interval on LSA
		{			
			string left_pattern = pattern.substr(start_pos, anchor_pos-start_pos+1);
			reverse(left_pattern.begin(), left_pattern.end());				
			
			INT rev_pair1=0,rev_pair2=0;
			rev_pattern_matching ( left_pattern, whole_string, LSA, LLCP, lrmq, g, rev_pair1,rev_pair2 );  
	
			if(rev_pair1 > rev_pair2)
			{	
				start_pos = start_pos + 1; 
				continue;
			}
			
			INT anchor_pos_to_return=anchor_pos;
			for(INT i = rev_pair1; i <= rev_pair2; i++ ) //this can be an interval of size 10^10 and only one occurrence is valid.
			{
				INT index = whole_string_len-1-LSA[i];	
				INT index_to_return=index;
				anchor_pos=anchor_pos_to_return;
				
					index++; anchor_pos++;
				
					while ( ( anchor_pos < start_pos+ell ) && ( index < whole_string_len ) && ( whole_string[index] == pattern[anchor_pos] ) )
					{
						anchor_pos++;						
						index++;
					
					}
					
					if(anchor_pos==start_pos+ell)
					{
						pair<INT,INT> my_pair(anchor_pos_to_return,invRSA[whole_string_len-1-LSA[i]]);			
						pattern_answers.push_back(my_pair);
					}
				
				
			}
		}				
	
     start_pos = start_pos + 1; 
  }
  
  delete []f;
}


/* Constructs the right compacted trie given the anchors and the SA of the whole string */
void right_compacted_trie ( unordered_map<INT,INT> &anchors, INT * SA, INT * LCP, INT n, INT * RSA, INT * RLCP, unordered_map<INT,INT> &invRSA, INT g )
{
	INT ii = 0; 
	INT minLCP = n;

	for( INT i = 0; i < n; i++ ) 
	{
		
		auto it = anchors.find( SA[i] );
		if( it != anchors.end() )
		{
			RSA[ii] = SA[i];		
			invRSA[RSA[ ii ]] = ii;
			if ( ii == 0 )	RLCP[ii] = 0; 
			else
			{
				if ( SA[i-1] == RSA[ii-1] ) 
					RLCP[ii] = LCP[i];	
				else
					RLCP[ii] = std::min(minLCP, LCP[i]);
			}
  			
			minLCP = n; 
			ii++;
		}
		else 
		{
			if ( LCP[i] < minLCP )
				minLCP = LCP[i];
		}
	}
}

/* Constructs the left compacted trie given the anchors and the SA of the whole string */
void left_compacted_trie ( unordered_map<INT,INT> &anchors, INT * SA, INT * LCP, INT n, INT * RSA, INT * RLCP, unordered_map<INT,INT> &invRSA, INT g )
{
	INT ii = 0;
	INT minLCP = n;
	for( INT i = 0; i < n; i++ ) 
	{
		
		auto it = anchors.find( ( n - 1 ) - SA[i] );
		if( it != anchors.end() )
		{
			RSA[ii] = SA[i];		
			invRSA[(n-1) - RSA[ii]] = ii;	
			if ( ii == 0 )	RLCP[ii] = 0; 	
			else
			{
				if ( SA[i-1] == RSA[ii-1] ) 
					RLCP[ii] = LCP[i];	
				else
					RLCP[ii] = std::min(minLCP, LCP[i]);	
			}
  			
			
			minLCP = n; 
			ii++;
		}
		else 
		{
			if ( LCP[i] < minLCP )
				minLCP = LCP[i];
		}
	}
}

/* Computes the length of lcs of two suffixes of two strings */
INT lcs ( string & x, INT M, string & y, INT l )
{
	if ( M < 0 ) return 0;
	INT yy = y.size();
	if ( l >= yy ) return 0;

	INT i = 0;
	while ( ( M - i >= 0 ) && ( l + i < yy ) )
	{
		if ( x[M-i] != y[l+i] )	break;
		i++;
	}
	return i;
}

/* Searching a list of strings using LCP from "Algorithms on Strings" by Crochemore et al. Algorithm takes O(m + log n), where n is the list size and m the length of pattern */
void rev_pattern_matching ( string & w, string & a, INT * SA, INT * LCP, std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> &rmq, INT n, INT &pair1, INT &pair2 )
{
	INT m = w.size(); 
	INT N = a.size(); 
	INT d = -1;
	INT ld = 0;
	INT f = n;
	INT lf = 0;
		
	while ( d + 1 < f )
	{
		INT i = (d + f)/2;
		INT revSA = N - 1 - SA[i];
		std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >>::iterator it;

		/* lcp(i,f) */
		INT lcpif;
		it = rmq.find(make_pair(i+1, f));
		lcpif = it -> second;
		
		/* lcp(d,i) */
		INT lcpdi;
		it = rmq.find(make_pair(d+1, i));
		lcpdi = it -> second;

		if ( ( ld <= lcpif ) && ( lcpif < lf ) )
		{
			d = i;
			ld = lcpif;
		}
		else if ( ( ld <= lf ) && ( lf < lcpif ) ) 	f = i;
		else if ( ( lf <= lcpdi ) && ( lcpdi < ld ) )
		{
			f = i;
			lf = lcpdi;
		}	
		else if ( ( lf <= ld ) && ( ld < lcpdi ) )	d = i;
		else
		{
			INT l = std::max (ld, lf);
			l = l + lcs ( a, revSA - l, w, l );
			if ( l == m ) //lower bound is found, let's find the upper bound
		        {
				INT e = i;
				while ( d + 1 < e )
				{
					INT j = (d + e)/2;
					
					/* lcp(j,e) */
					INT lcpje;
					it = rmq.find(make_pair(j+1, e));
					lcpje = it -> second;

					if ( lcpje < m ) 	d = j;
					else 			e = j;
				}

				/* lcp(d,e) */
				INT lcpde;
				it = rmq.find(make_pair(d+1, e));
				lcpde = it -> second;

				if ( lcpde >= m )	d = std::max (d-1,( INT ) -1 );

				e = i;
				while ( e + 1 < f )
				{
					INT j = (e + f)/2;
					
					/* lcp(e,j) */
					INT lcpej;
					it = rmq.find(make_pair(e+1, j));
					lcpej = it -> second;

					if ( lcpej < m ) 	f = j;
					else 			e = j;
				}

				/* lcp(e,f) */
				INT lcpef;
				it = rmq.find(make_pair(e+1, f));
				lcpef = it -> second;

				if ( lcpef >= m )	f = std::min (f+1,n);

				pair1=d+1;
				pair2=f-1;
				return;
			}
			else if ( ( l == N - SA[i] ) || ( ( revSA - l >= 0 ) && ( l != m ) && ( a[revSA - l] < w[l] ) ) )
			{
				d = i;
				ld = l;
			}
			else
			{
				f = i;
				lf = l;
			}	

		}
	}
	
	pair1=d+1;
	pair2=f-1;
	return;
}



/* A fast RMQ data structure for Manber & Myers algorithm: construct the recursion tree over array arr; at each node (i,j) of the tree we maintain the minimum value in arr[i+1..j] */
INT fast_RMQ ( INT * arr, INT L, INT R, INT n, std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> &rmq )
{
	if ( R - L > 1 ) // Internal nodes
	{
		INT M = (L+R)/2;
		INT a = fast_RMQ ( arr, L, M, n, rmq ); //Recurse on the left
		INT b = fast_RMQ ( arr, M, R, n, rmq ); //Recurse on the right
		INT value = std::min(a, b); 	     //This is the minimum of arr[L+1..M]

		pair<INT,INT> p(L+1, R);
		pair<pair<INT,INT>,INT> np(p, value);
		rmq.insert(np);
	
		return value;
	}
	else  // Base case: lead nodes
	{
		INT value;
		if ( R >= n ) value = 0;
		else value = std::min (arr[L+1], arr[R]);
		
		pair<INT,INT> p(L+1, R);
		pair<pair<INT,INT>,INT> np(p, value);
		rmq.insert(np);
		return value;
	}
}

/* Computes the length of lcp of two suffixes of two strings */
INT lcp ( string & x, INT M, string & y, INT l )
{
	INT xx = x.size();
	if ( M >= xx ) return 0;
	INT yy = y.size();
	if ( l >= yy ) return 0;

	INT i = 0;
	while ( ( M + i < xx ) && ( l + i < yy ) )
	{
		if ( x[M+i] != y[l+i] )	break;
		i++;
	}
	return i;
}

/* Searching a list of strings using LCP from "Algorithms on Strings" by Crochemore et al. Algorithm takes O(m + log n), where n is the list size and m the length of pattern */
void pattern_matching ( string & w, string & a, INT * SA, INT * LCP, std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> &rmq, INT n, INT &pair1, INT &pair2 )
{
	INT m = w.size(); //length of pattern
	INT N = a.size(); //length of string 
	INT d = -1;
	INT ld = 0;
	INT f = n;
	INT lf = 0;
		
	
	while ( d + 1 < f )
	{
		INT i = (d + f)/2;
		std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >>::iterator it;

		/* lcp(i,f) */
		INT lcpif;
		it = rmq.find(make_pair(i+1, f));
		if ( it == rmq.end() )
		{
			cout<<"1. RMQ failure! Not able to find"<<i+1<<","<<f<<endl;
			break;
		}
		else lcpif = it -> second;
		
		/* lcp(d,i) */
		INT lcpdi;
		it = rmq.find(make_pair(d+1, i));
		if ( it == rmq.end() )
		{
			cout<<"2. RMQ failure! Not able to find"<<d+1<<","<<i<<endl;
			break;
		}
		else lcpdi = it -> second;

		if ( ( ld <= lcpif ) && ( lcpif < lf ) )
		{
			d = i;
			ld = lcpif;
		}
		else if ( ( ld <= lf ) && ( lf < lcpif ) ) 	f = i;
		else if ( ( lf <= lcpdi ) && ( lcpdi < ld ) )
		{
			f = i;
			lf = lcpdi;
		}	
		else if ( ( lf <= ld ) && ( ld < lcpdi ) )	d = i;
		else
		{
			INT l = std::max (ld, lf);
			l = l + lcp ( a, SA[i] + l, w, l );
			if ( l == m ) //lower bound is found, let's find the upper bound
		        {
				
				INT e = i;
				while ( d + 1 < e )
				{
					INT j = (d + e)/2;
					
					/* lcp(j,e) */
					INT lcpje;
					it = rmq.find(make_pair(j+1, e));
					if ( it == rmq.end() )
					{
						cout<<"3. RMQ failure! Not able to find"<<j+1<<","<<e<<endl;
						break;
					}
					else lcpje = it -> second;

					if ( lcpje < m ) 	d = j;
					else 			e = j;
				}

				/* lcp(d,e) */
				INT lcpde;
				it = rmq.find(make_pair(d+1, e));
				if ( it == rmq.end() )
				{
					cout<<"4. RMQ failure! Not able to find"<<d+1<<","<<e<<endl;
					break;
				}
				else lcpde = it -> second;

				if ( lcpde >= m )	d = std::max (d-1,( INT ) -1 );

				e = i;
				while ( e + 1 < f )
				{
					INT j = (e + f)/2;
					
					/* lcp(e,j) */
					INT lcpej;
					it = rmq.find(make_pair(e+1, j));
					if ( it == rmq.end() )
					{
						cout<<"5. RMQ failure! Not able to find"<<e+1<<","<<j<<endl;
						break;
					}
					else lcpej = it -> second;

					if ( lcpej < m ) 	f = j;
					else 			e = j;
				}

				/* lcp(e,f) */
				INT lcpef;
				it = rmq.find(make_pair(e+1, f));
				if ( it == rmq.end() )
				{
					cout<<"6. RMQ failure! Not able to find"<<e+1<<","<<f<<endl;
					break;
				}
				else lcpef = it -> second;

				if ( lcpef >= m )	f = std::min (f+1,n);

				pair1 = d + 1;
				pair2 = f - 1;
				return;
				
				
			}
			else if ( ( l == N - SA[i] ) || ( ( SA[i]+l < N ) && ( l != m ) && ( a[SA[i]+l] < w[l] ) ) )
			{
				d = i;
				ld = l;
			}
			else
			{
				f = i;
				lf = l;
			}	

		}
	}

	pair1 = d + 1;
	pair2 = f - 1;
	return;
}


//ADAPTED from https://cp-algorithms.com/sequences/longest_increasing_subsequence.html#toc-tgt-6
vector<pair<INT,INT> > lis(vector<pair<INT,INT> > const& a) {
    INT n = a.size();
    vector<INT> d(n, 1), p(n, -1);
    for (INT i = 0; i < n; i++) {
        for (INT j = 0; j < i; j++) {
            if (a[j].second < a[i].second && d[i] < d[j] + 1) {
                d[i] = d[j] + 1;
                p[i] = j;
            }
        }
    }

    INT ans = d[0], pos = 0;
    for (INT i = 1; i < n; i++) {
        if (d[i] > ans) {
            ans = d[i];
            pos = i;
        }
    }

    vector<pair<INT,INT> > subseq;
    while (pos != -1) {
        subseq.push_back(a[pos]);
        pos = p[pos];
    }
    reverse(subseq.begin(), subseq.end());
    return subseq;
}

struct IndexValueCompare{
    inline bool operator() (const pair<INT, INT> &one, const pair<INT, INT> &another){
        return one.second < another.second;
    }
};

vector<pair<INT,INT> > lis_new(const vector<pair<INT, INT> > &sequence){
    vector<INT> parent(sequence.size());
    set<pair<INT, INT>, IndexValueCompare> s;
    for(INT i = 0; i < sequence.size(); ++i){
        pair<INT, INT> iv(i, sequence[i].second);
        if(i == 0){
            s.insert(iv);
            continue;
        }
        auto index = s.lower_bound(iv);
        if(index != s.end()){
            if(sequence[i].second < sequence[index->first].second){
                if(index != s.begin()) {
                    parent[i] = (--index)->first;
                    index++;
                }
                s.erase(index);
            }
        } else{
            parent[i] = s.rbegin()->first;
        }
        s.insert(iv);
    }
    vector<pair<INT,INT> > result(s.size());
    int index = s.rbegin()->first;
    for(auto iter = s.rbegin(); iter != s.rend(); index = parent[index], ++iter){
        result[distance(iter, s.rend()) - 1] = sequence[index];
    }
    return result;
}

  
//adapted from https://www.geeksforgeeks.org/how-to-find-index-of-a-given-element-in-a-vector-in-cpp/
void getIndex(vector<INT> &v, INT K, vector<double> &indices)
{
    auto it = find(v.begin(), v.end(), K);
 
    // If element was found
    if (it != v.end()) 
    {
     
        // calculating the index
        // of K
        int index = it - v.begin();
  
		indices.push_back((double)index);
    }
    else {
        // If the element is not
        // present in the vector
 
		indices.push_back(-1.0);
    }
}


int main(int argc, char** argv)
{
 
 if(argc!=7)
 {
        cout<<"Wrong arguments\n";
 	cout<<"./BDA_Search_v2 data_file ell pattern_file K tau delta \n";
 	exit(-1);
 }
 
 ifstream is;
 is.open (argv[1], ios::in | ios::binary);
 
 INT ell;
 std::stringstream(argv[2])>>ell;
 
 
 ifstream is2; 
 
 is2.open (argv[3], ios::binary);
 
 INT K;
 stringstream(argv[4])>>K;
 
 INT tau;
 stringstream(argv[5])>>tau;

 INT delta;
 stringstream(argv[6])>>delta; 
 
 cout<<"Parameters: ell = "<<ell<<", K="<<K<<", tau="<<0<<" delta="<<0<<endl;
 cout<<"(In future versions there will be support for tau>0 and delta>0)\n";
 
  vector<pair<INT, vector<unsigned char> > > all_strings;
  vector<unsigned char> input_seq_vec;
  INT id=0;
  char c=0; 
  while (is.read(reinterpret_cast<char*>(&c), 1))     
  {
      if(c=='\n')
      {
		if(input_seq_vec.empty())
			break; 
		vector<unsigned char> v(input_seq_vec);
		pair<INT, vector<unsigned char> > p(id, v);
		all_strings.push_back(p);
		input_seq_vec.clear();
		id++;
      }
      else
      {
        input_seq_vec.push_back((unsigned char)c);
      }
  }
  is.close();      
  
   
  INT *num_anchors_per_string=(INT *)malloc(id*sizeof(INT));
  //key is the anchor plus offset and value string_id 
  // to each anchor 
 unordered_map<INT, INT> global_anchors;
  INT total_length_of_previous_strings=0;
  
  //for each string in the database
  for(INT j=0;j<id;++j)
  {
	  
	
	  unordered_set<INT> anchors;
	  string current_string(all_strings[j].second.begin(),all_strings[j].second.end()); 
	  	
  
	  
      fast_anchors(current_string, anchors, ell);
	
	  for(auto & it : anchors)
	  {
		  INT anchor = it + total_length_of_previous_strings; 
		  //we add anchor and string_id
		  global_anchors.insert(make_pair(anchor,j));
		   
	  }
	
	  //how many anchors each string has
	  num_anchors_per_string[j]=anchors.size();
	  
	  total_length_of_previous_strings +=current_string.size();		
  }
  
  
  unsigned char *my_whole_string=(unsigned char*) malloc(total_length_of_previous_strings*sizeof(INT));
  
  for(INT i=0;i<total_length_of_previous_strings;)
  {
	for(INT j=0;j<id;++j)  
	{
		string current_string(all_strings[j].second.begin(),all_strings[j].second.end()); 
		for(INT k=0;k<current_string.size();++k)
			my_whole_string[i++]=current_string[k];  
	}
  }
  string whole_string( reinterpret_cast<char const*>(my_whole_string), total_length_of_previous_strings ) ;

  free(my_whole_string);
  
  
  INT g=global_anchors.size();
  

  INT * SA;
  INT * LCP;
  INT * invSA;

  INT n = whole_string.size();
  unsigned char * seq = ( unsigned char * ) whole_string.c_str();

  /* Compute the suffix array */
  SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
  if( ( SA == NULL) )
  {
  	fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
        return ( 0 );
  }
	
	#ifdef _USE_64
  	if( divsufsort64( seq, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
	#endif

	#ifdef _USE_32
  	if( divsufsort( seq, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
	#endif
	
  
  /*Compute the inverse SA array */
  invSA = ( INT * ) calloc( n , sizeof( INT ) );
  if( ( invSA == NULL) )
  {
  	fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
        return ( 0 );
  }
  for ( INT i = 0; i < n; i ++ )
  {
  	invSA [SA[i]] = i;
  }

  /* Compute the LCP array */
  LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
  if( ( LCP == NULL) )
  {
  	fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
        return ( 0 );
  }
  if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
  {
	  fprintf(stderr, " Error: LCP computation failed.\n" );
          exit( EXIT_FAILURE );
  }
  
 
  /* Constructing right and left compacted tries */
  INT * RSA;
  INT * RLCP;
  unordered_map<INT,INT> invRSA;

  RSA = ( INT * ) malloc( ( g ) * sizeof( INT ) );
  if( ( RSA == NULL) )
  {
  	fprintf(stderr, " Error: Cannot allocate memory for RSA.\n" );
        return ( 0 );
  }

  RLCP = ( INT * ) malloc( ( g ) * sizeof( INT ) );
  if( ( RLCP == NULL) )
  {
  	fprintf(stderr, " Error: Cannot allocate memory for RLCP.\n" );
        return ( 0 );
  }

 
  right_compacted_trie ( global_anchors, SA, LCP, n, RSA, RLCP, invRSA, g );  
  
  //To get rid of global_anchors
  //I keep ids. ids[p] = string_id
  INT *ids=(INT *)malloc((g)*sizeof(INT));
  
 
  for(INT i=0;i<g;++i)
  {

	  ids[i]=global_anchors[RSA[i]];
  }
  
  /* We reverse the string for the left direction and also overwrite all other DSs */
  reverse(whole_string.begin(), whole_string.end());
  
  
  #ifdef _USE_64
  	if( divsufsort64( seq, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
	#endif
	#ifdef _USE_32
  	if( divsufsort( seq, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
	#endif
	
  
  for ( INT i = 0; i < n; i ++ )
  {
  	invSA [SA[i]] = i;
  }
  
  if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
  {
	  fprintf(stderr, " Error: LCP computation failed.\n" );
          exit( EXIT_FAILURE );
  }

  INT * LSA;
  INT * LLCP;
  unordered_map<INT,INT> invLSA;

  LSA = ( INT * ) malloc( ( g ) * sizeof( INT ) );
  if( ( LSA == NULL) )
  {
  	fprintf(stderr, " Error: Cannot allocate memory for LSA.\n" );
        return ( 0 );
  }
  LLCP = ( INT * ) malloc( ( g ) * sizeof( INT ) );
  if( ( LLCP == NULL) )
  {
  	fprintf(stderr, " Error: Cannot allocate memory for LLCP.\n" );
        return ( 0 );
  }

  left_compacted_trie ( global_anchors, SA, LCP, n, LSA, LLCP, invLSA, g );


  global_anchors.clear();
  
  /* After constructing the tries these DSs over the whole string are not needed anymore, our data structure must be of size O(g) */
  
  free ( SA );
  free ( LCP );
  free ( invSA );
  
  //After this whole_string will be ok (left to right)
  reverse(whole_string.begin(), whole_string.end());
 
    
  /* The following RMQ data structures are used for spelling pattern over the LSA and RSA */
  std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> lrmq;
  std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> rrmq;
  

  fast_RMQ ( LLCP, -1, g, g, lrmq ); 
  fast_RMQ ( RLCP, -1, g, g, rrmq ); 

	
  #if 0
  vector < INT > v ( g, 0 );
  for ( INT i = 0; i < g; i++ )
  {
	  v[i] = LLCP[i];
	  
  }
  rmq_succinct_sct<> lrmq(&v); //using sdsl
	
  
  for ( INT i = 0; i < g; i++ )
  {
	  v[i] = RLCP[i];
	  
  }
  rmq_succinct_sct<> rrmq(&v); //using sdsl

  util::clear(v);
  #endif
   
   	
	
  	invLSA.clear();
		
   
  vector<vector<unsigned char> > all_patterns;
  
  vector<unsigned char> new_vec;
  unsigned char d=0;
  unsigned char d2=0;
  while (is2.read(reinterpret_cast<char*>(&d), sizeof(d)))     
  {
	  if (!is2)
		cerr << "Read failed." << endl;
	
      if(d=='\n')
      {
		if(new_vec.empty())
			break;
		vector<unsigned char> v;
		for(auto itv=new_vec.begin();itv!=new_vec.end();++itv)
		{
			v.push_back(*itv);
		}		
		all_patterns.push_back(v);
		new_vec.clear();		
      }
      else
      {
        new_vec.push_back(d);		
      }
  }
  is2.close();      
  
 auto total_duration=0,total_duration2=0;
 
 INT pruned_with_tau=0;  
 
 INT cluster_id=0; 
 for(auto &it_pat : all_patterns)   
 {
	 /*                       QUERYING                               */
   
   string pattern(it_pat.begin(),it_pat.end());
      
   //keeps the ids of the top-K results w.r.t. edit distance 	
   list<INT> ground_truth_ids;	
   
   
   for(INT i=0;i<K;++i)
   {
	   ground_truth_ids.push_back(cluster_id*K+i);
	
   }
   
   cluster_id++;   
   
	/*  Estimate identity matching heuristic */
	
	unordered_map<INT, INT> anchors_and_window_ends_in_pattern;
	vector<pair<INT,INT> > pattern_answers2;

	//Consruct the input to LIS for each string - key is string_ID 
 	unordered_map<INT, vector<pair<INT,INT> > > all_LIS2;

    priority_queue<pair<double, INT> > pq_est;	
	
	   
	auto start_time = chrono::steady_clock::now();
	
	fast_anchors_pattern_end(whole_string, pattern, ell, LSA, LLCP, lrmq, RSA, RLCP, rrmq, g, pattern_answers2,anchors_and_window_ends_in_pattern, invRSA);     
 
    total_duration2 += chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-start_time).count();	
	
  
  for(auto &it : pattern_answers2)
  {
	  
	   
	  INT string_id;
	  string_id= ids[it.second]; 
	 
	  
        auto it3=all_LIS2.find(string_id);
	  
	    if(it3!=all_LIS2.end())
	    {
		
		  all_LIS2[string_id].push_back(make_pair(it.first, RSA[it.second]));
		  		
	    }		  
	    else
	    {
		
		  vector<pair<INT,INT> > v;
		  v.push_back(make_pair(it.first,RSA[it.second]));
		  		
		  all_LIS2.insert(make_pair(string_id,v));		  		  
	    }	  
  }
 
   
   
   for(auto &it : all_LIS2)
   {				 	  
		vector<pair<INT,INT> > ret=lis_new(it.second); 		
		
		//If LIS length < tau continue
		//Commented out as we now support only tau=0
		//if(ret.size()<tau)       
		//	continue;							


		vector<INT> e_is;
				
		
		for(auto &it2 : ret)
		{
			
			e_is.push_back(anchors_and_window_ends_in_pattern[it2.first]);
		}
		
		
		INT score=ell;
		
	   for(vector<INT>::iterator it2=e_is.begin();it2!=e_is.end();++it2)
		{
		
			auto it_next=it2;
			it_next++;
			if(it_next!=e_is.end())
			{
				
				if(*it2 + ell <=*it_next)
				{
					score+=ell;
		
				}
				else
				{
					score+=*it_next - *it2;
		
				}
			}
		}
	   
	   string cur_string(all_strings[it.first].second.begin(),all_strings[it.first].second.end());
	   
		
		double score_prod=-1.0*score;
		
		if(pq_est.size()<K) //if I have seen <K elements
		{
			pq_est.push(make_pair(score_prod,it.first));
		
		}
		else
		{
			//if the pq top (which is the minimum) is larger than the current minimum , remove it and add the current minimum 
			if(pq_est.top().first>(score_prod))
			{		
				pq_est.pop();
				pq_est.push(make_pair(score_prod,it.first));
			}			
		}		
   }
	total_duration += chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-start_time).count();	
	
	
	
   list<INT> est_ids;
      
   vector<INT> edit_distances_of_answers;   
   while(!pq_est.empty())
	{
		pair<double,INT> p=pq_est.top();
		
		pair<INT, vector<unsigned char> >  v=all_strings[p.second];
		string cur_string(v.second.begin(),v.second.end());
		
		INT ed=EditDistDP(pattern,cur_string);
		edit_distances_of_answers.push_back(ed);
		
		est_ids.push_front(p.second);
		
		cout<<"Edit distance: "<<ed<<"\n";
		cout<<"String       : "<<cur_string<<"\n";
		
		
		pq_est.pop();
	}
	
	if(est_ids.size()<K)
	{
		cout<<"Try a smaller K.\n";
		exit(-1);
	}
	    
	
 } 
 
 cout << "Elapsed total time (ms) : "<<total_duration<<endl;
 cout<< "Avg time per query (ms) : "<<(double)total_duration/all_patterns.size()<<endl;
   

 	free ( RSA );
  	free ( RLCP );
  	free ( LSA );
  	free ( LLCP );
  	free(ids);
	free(num_anchors_per_string);		
	invRSA.clear();
	    
	return 0;
}
