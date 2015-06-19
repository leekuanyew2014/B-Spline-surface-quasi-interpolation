#include<windows.h>
#include<stdio.h>
#include<vector>
#include<iostream>
#include<stdlib.h>
#include<fstream>
#include<math.h>

#define DRAW_BSPLINE_STEP 10
#define E 0.000000000000000000000001
#define MAX 20 

using namespace std;

struct Point3D
{
	double x;
	double y;
	double z;
};

class Point3
{
public:
	double x;
	double y;
	double z;
	int repeat;

	Point3(){};
	Point3(double xx, double yy, double zz) : x(xx), y(yy), z(zz), repeat(1) {};

//	Point(double xx, double yy, int re) : x(xx),y(yy),repeat(re) {}
	Point3(const Point3 &p, int re) : x(p.x), y(p.y), z(p.z), repeat(re) {};

	Point3 operator * (double number)
	{
		return Point3(x * number, y * number, z * number); 
	};

	bool operator == (Point3 p1)
	{
		if ((p1.x == x) && (p1.y == y) && (p1.z == z))
		{
			return TRUE;
		}
		else
			return FALSE;
	};

	bool operator != (Point3 p1)
	{
		if ((p1.x == x) && (p1.y == y) && (p1.z == z))
		{
			return FALSE;
		}
		else
			return TRUE;
	};

	int getRepeat()
	{	
		return this->repeat;
	};

	void setRepeat(int re)
	{
		this->repeat = re;
	};	
};

typedef struct {
	int		num_pnt;	// the number of points
    Point3D *Point;		// the point array
	int		closed;		// 1=closed 0=open
} Polygon3D;

vector<vector<Point3> > ShapePoints3D;//型值点

vector<vector<Point3> > ControlPoints3D;//控制顶点

vector<Point3> ShapePointsWithRepeat3;//每条曲线上的型值点

vector<Point3> AllPoints3;//全部的点

vector<double> NodeVector3;//节点向量

vector<Point3> ControlPointsQuasi3;//每条曲线的控制顶点

vector<Point3> BSplinePointsQuasi3;//曲线上的点

vector<vector<Point3> > SurfaceControlPoints(vector<vector<Point3> > ShapePoints3D);//得到曲面的控制顶点

double GetElement3D(int num, double t, vector<Point3> ControlPointsQuasi3, vector<double> NodeVector3);//得到值

vector<double> GetNodeVectorUnUniform3D(vector<Point3> ShapePointsWithRepeat3);//得到节点向量

Point3 GetdeBoorValue3D(int h, double t, vector<Point3> ControlPointsQuasi3);//根据de-boor算法得到函数值

vector<Point3> GetControlPnts3D(vector<Point3> ShapePointsWithRepeat3, vector<double> NodeVector3);//计算控制顶点

vector<Point3> GetBSplinePntsQuasi3D();//得到曲线上的点，在曲面上已经淘汰了

/**
 * 计算矩阵src的模
 */ 
double calculate_A(double src[][MAX], int n ) 
{ 
    int i,j,k,x,y; 
    double tmp[MAX][MAX], t; 
    double result = 0.0; 
       
    if( n == 1 ) 
    { 
        return src[0][0]; 
    } 
       
    for( i = 0; i < n; ++i ) 
    { 
        for( j = 0; j < n - 1; ++j ) 
        { 
            for( k = 0; k < n - 1; ++k ) 
            { 
                x = j + 1; 
                y = k >= i ? k + 1 : k; 
                   
                tmp[j][k] = src[x][y]; 
            } 
        } 
           
        t = calculate_A( tmp, n - 1 ); 
           
        if( i % 2 == 0 ) 
        { 
            result += src[0][i] * t; 
        } 
        else 
        { 
            result -= src[0][i] * t; 
        } 
    } 
   
    return result; 
} 
   
/**
 * 计算伴随矩阵
 */ 
void calculate_A_adjoint( double src[MAX][MAX], double dst[MAX][MAX], int n ) 
{ 
    int i, j, k, t, x, y; 
    double tmp[MAX][MAX]; 
   
    if( n == 1 ) 
    { 
        dst[0][0] = 1; 
        return; 
    } 
       
    for( i = 0; i < n; ++i ) 
    { 
        for( j = 0; j < n; ++j ) 
        { 
            for( k = 0; k < n - 1; ++k ) 
            { 
                for( t = 0; t < n - 1; ++t ) 
                { 
                    x = k >= i ? k + 1 : k ; 
                    y = t >= j ? t + 1 : t; 
                       
                    tmp[k][t] = src[x][y]; 
                } 
            } 
               
            dst[j][i]  =  calculate_A( tmp, n - 1 ); 
               
            if( ( i + j ) % 2 == 1 ) 
            { 
                dst[j][i] = -1*dst[j][i]; 
            } 
        } 
    } 
} 
   
/**
 * 得到逆矩阵
 */ 
int calculate_A_inverse( double src[MAX][MAX], double dst[MAX][MAX], int n ) 
{ 
    double A = calculate_A( src, n ); 
    double tmp[MAX][MAX]; 
    int i, j; 
   
    if ( fabs( A - 0 ) <= E ) 
    { 
        printf("不可能有逆矩阵！\n"); 
        return 0; 
    } 
   
    calculate_A_adjoint( src, tmp, n );   
   
    for( i = 0; i < n; ++i )   
    {   
        for( j = 0; j < n; ++j )   
        {   
            dst[i][j] = (double)( tmp[i][j] / A ); 
        }   
    } 
		   
    return 1; 
} 

vector<vector<Point3> > SurfaceControlPoints(vector<vector<Point3> > ShapePoints3D)
{
	ofstream fout("ControlPoints_stator.txt");
	ofstream f("ShapePoints_stator.txt");
	int i, j, l;
	vector<vector<Point3> > tmpPoints;//tmp数组，用来保存每条曲线的控制顶点
	vector<vector<Point3> > ControlPoints3D;//返回值，计算的是最终的控制顶点
	for (i = 0; i < ShapePoints3D.size(); i++)
	{
		vector<Point3> ShapePointsWithRepeat3 = ShapePoints3D[i];
		f<<ShapePointsWithRepeat3.size()<<endl;
		for (l = 0; l < ShapePointsWithRepeat3.size(); l++)
		{
			f<<ShapePointsWithRepeat3[l].x<<" "<<ShapePointsWithRepeat3[l].y<<" "<<ShapePointsWithRepeat3[l].z<<endl;
		}
		f<<endl;
		vector<double> NodeVector3 = GetNodeVectorUnUniform3D(ShapePointsWithRepeat3);//先得到节点向量
		vector<Point3> ControlPointsQuasi3 = GetControlPnts3D(ShapePointsWithRepeat3, NodeVector3);//得到每条曲线的控制顶点
		fout<<ControlPointsQuasi3.size()<<endl;
		for (j = 0; j < ControlPointsQuasi3.size(); j++)
		{
			fout<<ControlPointsQuasi3[j].x<<" "<<ControlPointsQuasi3[j].y<<" "<<ControlPointsQuasi3[j].z<<endl;
		}
		fout<<endl;
		tmpPoints.push_back(ControlPointsQuasi3);//将控制顶点保存
	}
	for (i = 0; i < tmpPoints[0].size(); i++)
	{
		vector<Point3> ShapePointsWithRepeat3;//记录每一列的控制顶点
		for (j = 0; j < tmpPoints.size(); j++)
		{
			ShapePointsWithRepeat3.push_back(tmpPoints[j][i]);
		}
		vector<double> NodeVector3 = GetNodeVectorUnUniform3D(ShapePointsWithRepeat3);
		vector<Point3> ControlPointsQuasi3 = GetControlPnts3D(ShapePointsWithRepeat3, NodeVector3);
		ControlPoints3D.push_back(ControlPointsQuasi3);
	}
	return ControlPoints3D;
}

double GetElement3D(int num, double t, vector<Point3> ControlPointsQuasi3, vector<double> NodeVector3)//num is BSpline index,t is canshu
{
	int k = 4;

	int size = ControlPointsQuasi3.size() - 1;
	
	if (fabs(t) < E)
	{
		if (num == 0)
		{
			return 1.0;
		}
		else 
		{	
			return 0.0;
		}
	}
	
	if (fabs(t - NodeVector3[size + 1]) < E)
	{
		if (num == size)
		{
			return 1.0;
		}
		else
		{
			return 0.0;
		}
	}
	

	if (t <= NodeVector3[num] || t >= NodeVector3[num + k])
	{
		return 0.0;
	}

	int id;

	for (int i = num; i < num + k; i++)
	{
		if (t >= NodeVector3[i] && t < NodeVector3[i + 1])
		{
			id = i;
			break;
		}
	}

	Point3* points3 = new Point3[k]; 
	
	int j;
	int index = 0;
	for (j = id - k + 1; j <= id; j++)
	{
		if (j == num)
		{
			points3[index].x = 1.0;
		}
		else
		{
			points3[index].x = 0.0;
		}
		index++;
	}
	int In = id + 1 - k;
	int l;
	double s;
	double den;
	for (l = 1; l < k; l++) //k次循环
	{
		for (j = k - 1; j >= l; j--)
		{
			den = NodeVector3[In + j + k - l] - NodeVector3[In + j];
			if (fabs(den) < E)
			{
				s = 0.0;
			}
			else
			{
				s = (t - NodeVector3[In + j]) / den;
			}
			points3[j].x = ((1.0 - s) * points3[j - 1].x + s * points3[j].x); 
		}
	}

	double returnPoint = points3[k - 1].x;
	delete[] points3;

	return returnPoint;
}

vector<double> GetNodeVectorUnUniform3D(vector<Point3> ShapePointsWithRepeat3)
{
	int num = ShapePointsWithRepeat3.size();//这里指的是控制顶点的数量减2
	
	vector<double> NodeVector3;
	
	if (num < 3)
	{
		return NodeVector3;
	}
	
	int size = num + 2;

	double* distance = new double[num - 1];
	double totalDistance = 0.0; //总弦长

	for (int i = 0; i < num - 1; i++)
	{
		double x1 = ShapePointsWithRepeat3[i].x;
		double x2 = ShapePointsWithRepeat3[i + 1].x;
		double y1 = ShapePointsWithRepeat3[i].y;
		double y2 = ShapePointsWithRepeat3[i + 1].y;
		double z1 = ShapePointsWithRepeat3[i].z;
		double z2 = ShapePointsWithRepeat3[i + 1].z;
		distance[i] = sqrt(pow(x1 - x2, 2)+pow(y1 - y2, 2) + pow(z1 - z2, 2));//这里本来有两个sqrt的
		totalDistance += distance[i];
	}

	for (i = 0; i < 4; i++)
	{
		NodeVector3.push_back(0.0);
	}

	double temp; 
	for (i = 4; i < size; i++)
	{
		temp = NodeVector3[i - 1] + (distance[i - 4] / totalDistance);
		NodeVector3.push_back(temp);
	}

	for (i = size; i < size + 4; i++)
	{
		NodeVector3.push_back(1.0);
	}
	return NodeVector3;
}

Point3 GetdeBoorValue3D(int h, double t, vector<Point3> ControlPointsQuasi3)//yangxunnian
{
	int k = 4; //曲线的阶数
	
	Point3* points3 = new Point3[k]; 

	int In = h + 1 - k;

	int j;
	int index = 0;
	for (j = h - k + 1; j <= h; j++)
	{
		points3[index].x = ControlPointsQuasi3[j].x;
		points3[index].y = ControlPointsQuasi3[j].y;
		points3[index].z = ControlPointsQuasi3[j].z;
		index++;
	}

	int l;
	double s;
	double den;
	for (l = 1; l < k; l++) //k次循环
	{
		for (j = k - 1; j >= l; j--)
		{
			den = NodeVector3[In + j + k - l] - NodeVector3[In + j];
			if (fabs(den) < E)
			{
				s = 0.0;
			}
			else
			{
				s = (t - NodeVector3[In + j]) / den;
			}
			points3[j].x = ((1.0 - s) * points3[j - 1].x + s * points3[j].x); 
			points3[j].y = ((1.0 - s) * points3[j - 1].y + s * points3[j].y);
			points3[j].z = ((1.0 - s) * points3[j - 1].z + s * points3[j].z);
		}
	}

	Point3 returnPoint = points3[k - 1];
	delete[] points3;

	return returnPoint;
}

vector<Point3> GetControlPnts3D(vector<Point3> ShapePointsWithRepeat3, vector<double> NodeVector3)
{
	vector<Point3> ControlPointsQuasi3;
	int num = ShapePointsWithRepeat3.size();

	if (num < 3)
	{	
		return ControlPointsQuasi3;
	}

	int i, j, k, l;

	Point3 point3(0.0, 0.0, 0.0);

	//int size = num + 2;//在这里多了两个控制顶点，是加上两个边界条件

	int size = num + 2;

	for (i = 0; i < size; i++)
	{
		ControlPointsQuasi3.push_back(point3);
	}	

	int n = 5;//fangchengzugeshu
	double matrixX[MAX][MAX];//coef matrix x

	double inverseMatrixX[MAX][MAX];//inverse matrix X

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			matrixX[i][j] = 0.0;
			inverseMatrixX[i][j] = 0.0;
		}
	}
	
	double *x = new double[n];
	/*
	x[0] = NodeVector[4 + 1];	
	x[2] = NodeVector[4 + 2]; 
	x[4] = NodeVector[4 + 3];
	x[6] = NodeVector[4 + 4];
	x[8] = NodeVector[4 + 5];

	x[1] = 1.0 / 2 * (x[0] + x[2]);
	x[3] = 1.0 / 2 * (x[2] + x[4]);
	x[5] = 1.0 / 2 * (x[4] + x[6]);
	x[7] = 1.0 / 2 * (x[6] + x[8]);
	
	for (j = 0; j < n; j++)
	{
		for (l = 0; l < n; l++)
		{
			matrixX[j][l] = GetElement5(l, x[j]);
		}
	}
				
	calculate_A_inverse(matrixX, inverseMatrixX, n);
	
	*/
	Point3 *points3 = new Point3[n];

	for (i = 2; i < size - 2; i++)
	{
		x[0] = NodeVector3[i + 1];	
		x[2] = NodeVector3[i + 2]; 
		x[4] = NodeVector3[i + 3];
		//x[6] = NodeVector[i + 4];
		//x[8] = NodeVector[i + 5];

		x[1] = 1.0 / 2 * (x[0] + x[2]);
		x[3] = 1.0 / 2 * (x[2] + x[4]);
		//x[5] = 1.0 / 2 * (x[4] + x[6]);
		//x[7] = 1.0 / 2 * (x[6] + x[8]);
		
		for (j = 0; j < n; j++)
		{	
			for (l = 0; l < n; l++)
			{
				matrixX[j][l] = GetElement3D(i - 2 + l, x[j], ControlPointsQuasi3, NodeVector3);//这里是不用改的
			}
		}

		calculate_A_inverse(matrixX, inverseMatrixX, n);
		/*
		for (int k  = 0; k < n; k++)
		{
			points[k] = fabssin((x[k] - 3) / (num - 3) * 4 - 2);
		}
		*/
		/*
		points[0] = ShapePointsWithRepeat[i - 2];
		points[2] = ShapePointsWithRepeat[i - 1];
		points[4] = ShapePointsWithRepeat[i];
		points[1].x = 1.0 / 2 * (points[0].x + points[2].x);
		points[1].y = 1.0 / 2 * (points[0].y + points[2].y);
		points[3].x = 1.0 / 2 * (points[2].x + points[4].x);
		points[3].y = 1.0 / 2 * (points[2].y + points[4].y);
		*/
		points3[0] = ShapePointsWithRepeat3[i - 2];
		points3[2] = ShapePointsWithRepeat3[i - 1];
		points3[4] = ShapePointsWithRepeat3[i];
		points3[1].x = 1.0 / 2 * (points3[0].x + points3[2].x);
		points3[1].y = 1.0 / 2 * (points3[0].y + points3[2].y);
		points3[1].z = 1.0 / 2 * (points3[0].z + points3[2].z);
		points3[3].x = 1.0 / 2 * (points3[2].x + points3[4].x);
		points3[3].y = 1.0 / 2 * (points3[2].y + points3[4].y);
		points3[3].z = 1.0 / 2 * (points3[2].z + points3[4].z);
		
		if (i == 2)
		{
			for (j = 0; j < n; j++)
			{
				ControlPointsQuasi3[0].x += inverseMatrixX[0][j] * points3[j].x;
				ControlPointsQuasi3[0].y += inverseMatrixX[0][j] * points3[j].y;
				ControlPointsQuasi3[0].z += inverseMatrixX[0][j] * points3[j].z;
				
				ControlPointsQuasi3[1].x += inverseMatrixX[1][j] * points3[j].x;
				ControlPointsQuasi3[1].y += inverseMatrixX[1][j] * points3[j].y;
				ControlPointsQuasi3[1].z += inverseMatrixX[1][j] * points3[j].z;
				//ControlPointsQuasi[2].x += inverseMatrixX[2][j] * points[j].x;
				//ControlPointsQuasi[2].y += inverseMatrixX[2][j] * points[j].y;
				//ControlPointsQuasi[3].x += inverseMatrixX[3][j] * points[j].x;
				//ControlPointsQuasi[3].y += inverseMatrixX[3][j] * points[j].y;
			}
			for (k = 0; k < n; k++)
			{
				for (l = 0; l < n; l++)
				{
					printf("%lf  ", inverseMatrixX[k][l]);
				}
				printf("\n");
			}
		}
		if (i == size - 3)
		{
			for (j = 0; j < n; j++)
			{	
				//ControlPointsQuasi[num - 4].x += inverseMatrixX[5][j] * points[j].x;
				//ControlPointsQuasi[num - 4].y += inverseMatrixX[5][j] * points[j].y;
				//ControlPointsQuasi[num - 3].x += inverseMatrixX[4][j] * points[j].x;
				//ControlPointsQuasi[num - 3].y += inverseMatrixX[4][j] * points[j].y;
				ControlPointsQuasi3[size - 2].x += inverseMatrixX[3][j] * points3[j].x;
				ControlPointsQuasi3[size - 2].y += inverseMatrixX[3][j] * points3[j].y;
				ControlPointsQuasi3[size - 2].z += inverseMatrixX[3][j] * points3[j].z;
				
				ControlPointsQuasi3[size - 1].x += inverseMatrixX[4][j] * points3[j].x;
				ControlPointsQuasi3[size - 1].y += inverseMatrixX[4][j] * points3[j].y;
				ControlPointsQuasi3[size - 1].z += inverseMatrixX[4][j] * points3[j].z;
			}
			for (k = 0; k < n; k++)
			{
				for (l = 0; l < n; l++)
				{
					printf("%lf  ", inverseMatrixX[k][l]);
				}
				printf("\n");
			}
		}
		
		for (j = 0; j < n; j++)
		{
			ControlPointsQuasi3[i].x += inverseMatrixX[2][j] * points3[j].x;
			ControlPointsQuasi3[i].y += inverseMatrixX[2][j] * points3[j].y;
			ControlPointsQuasi3[i].z += inverseMatrixX[2][j] * points3[j].z;
		}
	}
	for (i = 0; i < n; i++)
	{
		for (j =0; j < n;j++)
		{
			cout<<inverseMatrixX[i][j]<<"  ";
		}
		cout<<endl;
	}
	//ControlPointsQuasi[0] = ControlPointsQuasi[num - 2];
	//ControlPointsQuasi[1] = ControlPointsQuasi[num - 1];
	return ControlPointsQuasi3;
}

vector<Point3> GetBSplinePntsQuasi3D()
{	
	//ControlPointsQuasi.push_back(ControlPointsQuasi[0]);
	//ControlPointsQuasi.push_back(ControlPointsQuasi[1]);
	//ControlPointsQuasi.push_back(ControlPointsQuasi[2]);
	vector<Point3> BSplinePointsQuasi3;
	int conNum = ControlPointsQuasi3.size();
 	
 	Point3 tempPoint3;
	Point3 cpoint3;
 	for (int i = 3; i <= conNum-1; i++)
 	{
 		double step = (double(NodeVector3[i+1] - NodeVector3[i])) / DRAW_BSPLINE_STEP;
		if (fabs(step) < E)
		{
			continue;
		}
 		double initStep = NodeVector3[i];
 		
 		for (int j = 1; j < DRAW_BSPLINE_STEP; j++)
 		{
 				tempPoint3 = GetdeBoorValue3D(i, initStep, ControlPointsQuasi3);
				cpoint3.x = tempPoint3.x;
				cpoint3.y = tempPoint3.y;
				cpoint3.z = tempPoint3.z;
 				BSplinePointsQuasi3.push_back(cpoint3);
 				initStep += step;							
 		}

		tempPoint3 = GetdeBoorValue3D(i, NodeVector3[i + 1], ControlPointsQuasi3);
		cpoint3.x = tempPoint3.x;
		cpoint3.y = tempPoint3.y;
		cpoint3.z = tempPoint3.z;
 		BSplinePointsQuasi3.push_back(cpoint3);
 	}
	return BSplinePointsQuasi3;
};

Polygon3D* InputPolygon3D(char filename[])
{
	FILE   *fp;
	int    i, n, n0, cls=0;
	int tmp;
	int row;
	int col;
	Polygon3D *plyg;
	Point3D *P;

	if((fp = fopen(filename, "r")) == NULL) exit(1);
	fscanf (fp, "%d", &row);
	col = 71;
	n = row * 2 * (col + 1);
	
	n0 = row * 2 * col;
	if (n0>0)  {
		P = (Point3D *)malloc(n0 * sizeof (Point3D ));
	//	P = new Point3D[100];
	}
	else 
		return (Polygon3D *)NULL;
	int index = 0;
	for (i = 0; i < n; i++)
	{
		if (i % (col + 1) == 0)
		{
			fscanf (fp, "%d", &tmp);
			continue;
		}
	    fscanf (fp, "%lf%lf%lf",&(P[index].x), &(P[index].y), &(P[index].z));
		printf("P[%d]:%lf %lf %lf\n",i,P[index].x,P[index].y,P[index].z);
		index++;
	}

	fscanf (fp, "%d", &cls);
	fclose (fp);

	plyg = (Polygon3D *)malloc(sizeof (Polygon3D ));
	plyg->num_pnt=n0;
	plyg->Point=P;
	plyg->closed=(cls > 0)? 1 : 0;

	return plyg;
}

int main()
{
	char filename[12];
	printf("Please input the data  file:\n");
	scanf("%s", filename);
	Polygon3D *plyg;
	plyg = (Polygon3D *)InputPolygon3D(filename);
	Point3 point(0.0, 0.0, 0.0);
	vector<Point3> points;
	vector<vector<Point3> > ShapePoints3D;
	int i, j;
	for (i = 0; i < plyg->num_pnt; i++)
	{
		point.x = plyg->Point[i].x;
		point.y = plyg->Point[i].y;
		point.z = plyg->Point[i].z;
		points.push_back(point);
		if (i % 71 == 70)
		{
			ShapePoints3D.push_back(points);
			points.clear();
		}
		
	}
	vector<vector<Point3> > ControlPoints3D;
	ControlPoints3D = SurfaceControlPoints(ShapePoints3D);
	ofstream fout("stator.txt");
	fout<<ControlPoints3D.size()<<endl;
	for (i = 0; i < ControlPoints3D.size(); i++)
	{
		fout<<ControlPoints3D[0].size()<<endl;
		for (int j = 0; j < ControlPoints3D[i].size(); j++)
		{
			fout<<ControlPoints3D[i][j].x<<" "<<ControlPoints3D[i][j].y<<" "<<ControlPoints3D[i][j].z;
			fout<<endl;
		}
		fout<<endl;
	}
	cin.get();
	cin.get();
	return 0;
}