/**
@brief Find all set partitions with a given number of blocks.
  
https://codereview.stackexchange.com/questions/1526/finding-all-k-subset-partitions/1944#1944

For an implementation in Ruby, see:
https://gist.github.com/mikejholly/759ee1747b5ddb5ba052

Knuth's Algorithm U. The Art of Computer Programming, Volume 4,
Fascicle 3B.
**/

typedef struct Partitions Partitions;

struct Partitions {
    tipId_t u; // universe
    int sizeU; // number of elements in universe
    int nsub;  // number of subsets: 0 < nsub < sizeU
    int nway;  // number of partitions with n subsets

    // Matrix of dimension nway X nsub.
    // part[i][j] is the j'th subset within the i'th partition.
    tipId_t **part; 
};

Partitions *Partitions_new(int
