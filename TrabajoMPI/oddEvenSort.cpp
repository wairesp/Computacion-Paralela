
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>

using namespace  std;

/// MERGER SORT

int merge(double *ina, int lena, double *inb, int lenb, double *out) {
    int i,j;
    int outcount=0;

    for (i=0,j=0; i<lena; i++) {
        while ((inb[j] < ina[i]) && j < lenb) {
            out[outcount++] = inb[j++];
        }
        out[outcount++] = ina[i];
    }
    while (j<lenb)
        out[outcount++] = inb[j++];

    return 0;
}

int domy_merge(double *a, int start, int end, double *b) {
    if ((end - start) <= 1) return 0;

    int mid = (end+start)/2;
    domy_merge(a, start, mid, b);
    domy_merge(a, mid,   end, b);
    merge(&(a[start]), mid-start, &(a[mid]), end-mid, &(b[start]));
    for (int i=start; i<end; i++)
        a[i] = b[i];

    return 0;
}

int my_merge(int n, double *a) {
    double b[n];
    domy_merge(a, 0, n, b);
    return 0;
}

///// END MERGE


/// PRINT ITERACIONES
void printstat(int rank, int iter, char *txt, double *la, int n) {
    printf("[%d] %s iter %d: <", rank, txt, iter);
    for (int j=0; j<n-1; j++)
        printf("%6.3lf,",la[j]);
    printf("%6.3lf>\n", la[n-1]);
}

void MPI_Pairwise_Exchange(int localn, double *locala, int sendrank, int recvrank,
                           MPI_Comm comm) {

    /*
     * the sending rank just sends the data and waits for the results;
     * the receiving rank receives it, sorts the combined data, and returns
     * the correct half of the data.
     */
    int rank;
    double remote[localn];
    double all[2*localn];
    const int mergetag = 1;
    const int sortedtag = 2;

    MPI_Comm_rank(comm, &rank);
    if (rank == sendrank) {
        MPI_Send(locala, localn, MPI_DOUBLE, recvrank, mergetag, MPI_COMM_WORLD);
        MPI_Recv(locala, localn, MPI_DOUBLE, recvrank, sortedtag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(remote, localn, MPI_DOUBLE, sendrank, mergetag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        merge(locala, localn, remote, localn, all);

        int theirstart = 0, mystart = localn;
        if (sendrank > rank) {
            theirstart = localn;
            mystart = 0;
        }
        MPI_Send(&(all[theirstart]), localn, MPI_DOUBLE, sendrank, sortedtag, MPI_COMM_WORLD);
        for (int i=mystart; i<mystart+localn; i++)
            locala[i-mystart] = all[i];
    }
}

int MPI_OddEven_Sort(int n, double *a, int root, MPI_Comm comm)
{
    int rank, size, i;
    double *local_a;

// get rank and size of comm
    MPI_Comm_rank(comm, &rank); //&rank = address of rank
    MPI_Comm_size(comm, &size);

    local_a = (double *) calloc(n / size, sizeof(double));


// scatter the array a to local_a
    MPI_Scatter(a, n / size, MPI_DOUBLE, local_a, n / size, MPI_DOUBLE,
                root, comm);
// sort local_a
    my_merge(n / size, local_a);

//odd-even part
    for (i = 1; i <= size; i++) {
        char* txt="before";
        if ((i + rank) % 2 == 0) {  // means i and rank have same nature
            if (rank < size - 1) {
                MPI_Pairwise_Exchange(n / size, local_a, rank, rank + 1, comm);
            }
        } else if (rank > 0) {
            MPI_Pairwise_Exchange(n / size, local_a, rank - 1, rank, comm);
        }

    }

    printstat(rank, i-1, "after", local_a, n/size);

// gather local_a to a
    MPI_Gather(local_a, n / size, MPI_DOUBLE, a, n / size, MPI_DOUBLE,
               root, comm);

    if (rank == root)
        printstat(rank, i, " all done ", a, n);

    return MPI_SUCCESS;
}

int main(int argc,char ** argv) {

    srand(time(0));

    int n;
    cout<<"\ningrese n: ";
    cin>> n ;

    double a[n];

    for (int i=0; i<n; i++){
        a[i] = 1 + (rand() % 20);
        // cout<<a[i]<<" ";
    }


    MPI_Init(&argc,&argv);

    MPI_Barrier(MPI_COMM_WORLD);
    double start, finish, elapsed;
    start = MPI_Wtime();
    MPI_OddEven_Sort(n, a, 0, MPI_COMM_WORLD);
    finish = MPI_Wtime();

    elapsed = finish - start;
    printf("Elapsed time= %e seconds", elapsed);
    MPI_Finalize();

    return 0;
}


