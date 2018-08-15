#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <algorithm>
#include <vector>
#include <queue>
#include <iostream>
#include <string>
#include <map>
#include <list>
#include <stack>
#include <set>
#include <deque>

using namespace std;

#define OUTPUT freopen("myfile.txt","w",stdout);
#define INPUT freopen("myInput.in","r",stdin);
#define DEBUG(a) cout<<a<<endl;
#define PI acos(-1.0)
#define MAX 100005
#define MOD 1000000007
#define EPS 1e-9
#define BIGGER(a,b) (a>=b ? a : b)
#define SMALLER(a,b) (a<=b ? a : b)
#define getInt(a) scanf("%d",&a);
#define getLong(a) scanf("%lld",&a);
#define pb push_back
#define ppb pop_back
#define setBit(a,n) a|=1<<n;
#define resetBit(a,n) a&=~(1<<n);
#define checkBit(a,n) ((a&(1<<n))!=0)
#define toggleBit(a,n) a^=1<<n;

#define INF 1000000000
#define DBL_NEG_INF -1000000000.0

string seq1;
string seq2;
double match_score; 
double mismatch_penalty; 
double gap_start; 
double gap_extend; 
double memo[5005][5005][2];
int backtrack_memo[5005][5005][2];
char aligned_line_1[10010];
char aligned_line_2[10010];

double weight_value(int index1, int index2, int len1, int len2);
double dp(int index1, int index2, int continuing_gap);
double get_alignment();

double weight_value(int index1, int index2, int len1, int len2)
{
    // near the left terminal
    int x1_left = index1;
    int x2_left = index2;

    int x_left = SMALLER(x1_left, x2_left);

    // near the right terminal

    int x1_right = len1 - index1 -1 ;
    int x2_right = len2 - index2 -1 ;

    int x_right = SMALLER(x1_right, x2_right);

    // pick the smaller one
    double x = SMALLER(x1_left,x1_right) * 1.0;

    // normalize
    x = (2.0*x) / ((len1+len2)/2.0);  /*  len1 may be different from len2 so we added them
                                      *  and also multiply x by 2, so that the value becomes
                                      *  properly normalized.
                                      *  
                                      *  The final division by 2 is a bit complex to explain,
                                      *  as we mapped the corners to the middle for complicated convenience
                                      *  we needed to divide it by 2, so that the actual distance from the 
                                      *  center is obtained.
                                      *  
                                      *  Hopefully the manuscript will explain in a better way
                                      */

    // obtain weight from a modified gaussian function

    double miu = 0.0;
    double sigma = 0.4;

    double weight = exp(-((x-miu)*(x-miu)/(2.0*sigma*sigma)));
    
    return weight;
}

double dp(int index1, int index2, int continuing_gap)
{
    /*
     *  Notation :
     *                  'X' : Sentinal => 4
     *                  'D' : Diagonal => 1
     *                  'H' : Horizontal => 2
     *                  'V' : Vertical => 3
     */

    double maxx = DBL_NEG_INF;

    if(backtrack_memo[index1][index2][continuing_gap]!=0)
    {
        return memo[index1][index2][continuing_gap];
    }

    if((index1 == seq1.length()) && (index2 == seq2.length()))
    {
        backtrack_memo[seq1.length()][seq2.length()][continuing_gap] = 4;
        return 0.0;
    }

    else if ((index1 == seq1.length()) && (index2 != seq2.length()))
    {                                                                                   
        maxx = weight_value(index1, index2, seq1.length(), seq2.length()) * ( continuing_gap==1 ? gap_extend : gap_start ) + dp(index1, index2+1, 1);
        backtrack_memo[index1][index2][continuing_gap] = 3;
    }

    else if ((index1 != seq1.length()) and (index2 == seq2.length()))
    {
        maxx = weight_value(index1, index2, seq1.length(), seq2.length()) * ( continuing_gap==1 ? gap_extend : gap_start ) + dp(index1+1, index2, 1);
        backtrack_memo[index1][index2][continuing_gap] = 2;
    }

    else
    {
        // match 
        double score1 = DBL_NEG_INF;

        if(seq1[index1] == seq2[index2])
        {
            score1 = weight_value(index1, index2, seq1.length(), seq2.length()) * match_score + dp(index1+1, index2+1, 0);
        }
            

        // mismatch
        
        double score2 = DBL_NEG_INF;

        if(seq1[index1] != seq2[index2])
        {
            score2 = weight_value(index1, index2, seq1.length(), seq2.length()) * mismatch_penalty + dp(index1+1, index2+1, 0);
        }

        // gap

        // add gap in seq1
        double score3 = weight_value(index1, index2, seq1.length(), seq2.length()) * ( continuing_gap==1 ? gap_extend : gap_start ) + dp(index1, index2+1, 1);

        // add gap in seq2
        double score4 = weight_value(index1, index2, seq1.length(), seq2.length()) * ( continuing_gap == 1 ? gap_extend : gap_start ) + dp(index1+1, index2, 1);


        if( (score1-maxx) > EPS )
        {
            maxx = score1;
            backtrack_memo[index1][index2][continuing_gap] = 1;
        }

        if( (score2-maxx) > EPS )
        {
            maxx = score2;
            backtrack_memo[index1][index2][continuing_gap] = 1;
        }

        if( (score3-maxx) > EPS )
        {
            maxx = score3;
            backtrack_memo[index1][index2][continuing_gap] = 3;
        }

        if( (score4-maxx) > EPS )
        {
            maxx = score4;
            backtrack_memo[index1][index2][continuing_gap] = 2;
        }

    }

    memo[index1][index2][continuing_gap] = maxx ;

    return maxx;
}

double get_alignment()
{
    /*
     *  Notation :
     *                  'X' : Sentinal => 4
     *                  'D' : Diagonal => 1
     *                  'H' : Horizontal => 2
     *                  'V' : Vertical => 3
     */

    int line_ptr = 0;

    int index1 = 0;
    int index2 = 0;
    int gap_extend = 0;

    double norm_value=0.0;

    while(backtrack_memo[index1][index2][gap_extend] != 4)
    {
        if(backtrack_memo[index1][index2][gap_extend] == 1)
        {
            aligned_line_1[line_ptr] = seq1[index1];
            aligned_line_2[line_ptr] = seq2[index2];

            index1 += 1;
            index2 += 1;

            gap_extend = 0;
        }

        else if(backtrack_memo[index1][index2][gap_extend] == 2)
        {
            aligned_line_1[line_ptr] = seq1[index1];
            aligned_line_2[line_ptr] = '-';

            index1 += 1; 

            gap_extend = 1;
        }

        else if(backtrack_memo[index1][index2][gap_extend] == 3)
        {
            aligned_line_1[line_ptr] = '-';
            aligned_line_2[line_ptr] = seq2[index2];

            index2 += 1;

            gap_extend = 1;
        }

        line_ptr++;
    }

    aligned_line_1[line_ptr] = '\0';
    aligned_line_2[line_ptr] = '\0';

    for(int i=0;i<line_ptr;i++)
    {
        norm_value += weight_value(i, i, line_ptr, line_ptr);
    }

    if(fabs(norm_value)<EPS){
        norm_value = 1.0;
    }

    return norm_value;
        
}

int main(int argc, char** argv)
{
    seq1 = argv[1];
    seq2 = argv[2];
    
    match_score = atof(argv[3]); 
    mismatch_penalty = atof(argv[4]); 
    gap_start = atof(argv[5]); 
    gap_extend = atof(argv[6]); 

    int score_only = atoi(argv[7]); 

    
    memset(backtrack_memo,0,sizeof(backtrack_memo));

    double score = dp(0,0,0);
    double norm_value = get_alignment();

    if(score_only==1)
    {
        printf("%lf\n",score/norm_value);
    }

    
    else
    {        
        printf("%lf\n%s\n%s\n",score/norm_value,aligned_line_1,aligned_line_2);
    }
        
    

    return 0;
}