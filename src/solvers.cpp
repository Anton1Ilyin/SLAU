#include"matrix.cpp"

std::vector<double> YacMult(CSR_matrix mat,std::vector<double> vec)
{
    std::vector<double> res(vec.size());
    for (int i=0; i<vec.size(); i++)
        {
            for(int j=mat.atRow(i);j<mat.atRow(i+1);j++)
            {
                res[i]+=mat.atVal(j)*vec[mat.atCol(j)];
            }
        }
    return res;
}
CSR_matrix invdiag(CSR_matrix A)
{
    int n=A.RowSize()-1;
    std::vector<double> diag(n);
    std::vector<unsigned> col(n), row(n + 1);
    for(int i=0; i<n; i++)
    {
        for(int j=A.atRow(i); j<A.atRow(i+1); j++)
        {
            if(i==A.atCol(j))
            {
                diag[i]=1/A.atVal(j);
                col[i]=i;
                row[i+1]=i+1;
            }
        }
    } 
}

std::vector<double> GZ(CSR_matrix& A, std::vector<double>& b, std::vector<double> x0, double breakpoint)
{
    std::vector<double> x=x0;
    CSR_matrix d=invdiag(A);
    double r=norm2(A*x-b);
    double z;
    while(breakpoint<r)
    {
        for(int i=0; i<x.size(); i++)
        {
            z=0.0;
            for(int j=A.atRow(i); j<A.atRow(i+1); j++)
            {
                if(i!=A.atCol(j))
                {
                    z+=A.atVal(j)*x[A.atCol(j)];
                }
            }
            x[i]=d.atVal(i)*(b[i]-z);
        }
        r=norm2(A*x-b);
    }
    return x;
};

std::vector<double> Yacoby(CSR_matrix& A, std::vector<double>& b, std::vector<double> x0, double breakpoint)
{
    std::vector<double> x=x0;
    double r=norm2(A*x-b);
    while(breakpoint<r)
    {
        x=(b-YacMult(A,x));
        r=norm2(A*x-b);
    }
    return x;
};

std::vector<double> Ineration(CSR_matrix& A, std::vector<double>& b, double tau, std::vector<double> x0, double breakpoint)
{
    std::vector<double> x=x0;
    std::vector<double> r=b-A*x;
    while(breakpoint<norm2(r))
    {
        x=x+r*tau;
        r=b-A*x;
    }
    return x;
};


std::vector<double> Chebyshev(CSR_matrix& A, std::vector<double>& b, double lambda_min, double lambda_max, std::vector<double> x0, double breakpoint)
{
    std::vector<double> x=x0;
    std::vector<double> r=b-A*x;
    int Tn=128;
    int deg=7;
    std::vector<double> root(Tn);
    std::vector<double> index(Tn);
    index[0]=0;
    int k=Tn/4;
    index[Tn/2]=1;
    for(int i=2; i<deg; i++)
    {
        for(int j=0; j<Tn; j+=k)
        {
            index[j+k/2]=Tn/(k/2)-1-index[j];
        }
        k/=2;
    }
    for(int i=0; i<Tn; i++)
    {
        root[i]=(lambda_max+lambda_min)/2+(lambda_max-lambda_min)*cos((M_PI*(2*i+1))/(2*Tn))/2;
    }
    while(breakpoint<norm2(r))
    {
        for(int i=0; i<Tn; i++)
        {
            x=x+r/root[index[i]];
            r=b-A*x;
        }
    }
    return x;
};

std::vector<double> SGZ(CSR_matrix& A, std::vector<double>& b, std::vector<double> x0, double breakpoint)
{
    std::vector<double> x=x0;
    CSR_matrix d=invdiag(A);
    double r=norm2(A*x-b);
    double z;
    while(breakpoint<r)
    {
        for(int i=0; i<x.size(); i++)
        {
            z=0;
            for(int j=A.atRow(i); j<A.atRow(i+1); j++)
            {
                if(i!=A.atCol(j))
                {
                    z+=A.atVal(j)*x[A.atCol(j)];
                }
            }
            x[i]=d.atVal(i)*(b[i]-z);
        }
        for(int i=x.size()-1; i>0; i--)
        {
            z=0;
            for(int j=A.atRow(i); j<A.atRow(i+1); j++)
            {
                if(i!=A.atCol(j))
                {
                    z+=A.atVal(j)*x[A.atCol(j)];
                }
            }
            x[i]=d.atVal(i)*(b[i]-z);
        }
        r=norm2(A*x-b);
    }
    return x;
};

std::vector<double> Gradient(CSR_matrix& A, std::vector<double>& b, std::vector<double> x0, double breakpoint)
{
    std::vector<double> x=x0;
    std::vector<double> r=b-A*x;
    double tau=r*r/(r*(A*r));
    while(breakpoint<norm2(r))
    {
        x=x+r*tau;
        r=b-A*x;
        tau=r*r/(r*(A*r));
    }
    return x;
};