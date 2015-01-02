using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Media.Media3D;



namespace Wpf3DTube_NonStatic
{
    public class Plane
    {
        Vector3D prNormal;
        Point3D prPointOnPlane;
        PlaneEquation prPlaneEquation;
        public PlaneEquation Equation
        {
            get
            {
                return prPlaneEquation;
            }
            set
            {
                prPlaneEquation = value;
            }
        }
        public Vector3D Normal
        {
            get
            {
                return prNormal;
            }
            set
            {
                prNormal = value;
            }
        }
        public Point3D PointOnPlane
        {
            get
            {
                return prPointOnPlane;
            }
            set
            {
                prPointOnPlane = value;
            }
        }
        public static LineVector ReturnIntersectionLine_old(Plane A, Plane B)
        {
            //if normal is (a1,a2,a3) and point on plane is (x0,y0,z0)
            double x = 0;
            double y = 0;
            double d1 = A.Equation.Constantparameter;
            double d2 = B.Equation.Constantparameter;
            double b1 = A.Normal.Y;
            double b2 = B.Normal.Y;
            double a1 = A.Normal.X;
            double a2 = B.Normal.X;
            double c1 = A.Normal.Z;
            double c2 = B.Normal.Z;

            double denominator = a1 * b2 - a2 * b1;
            x = (d1 * b1 - d2 * b2) / (denominator);
            y = (a1 * d2 - a2 * d1) / (denominator);



            double vectorX = 0;
            double vectorY = 0;
            double vectorZ = 0;


            vectorX = c1 * b2 - c2 * b1;
            vectorY = a1 * c2 - a2 * c1;
            vectorZ = a2 * b1 - a1 * b2;

            Point3D aPoint = new Point3D(x, y, 0);
            Vector3D v = new Vector3D(vectorX, vectorY, vectorZ);
            return new LineVector(v, aPoint);
        }
        public static LineVector ReturnIntersectionLine(Plane A, Plane B)
        {
            //if normal is (a1,a2,a3) and point on plane is (x0,y0,z0)


            Point3D aPoint = SolveForOnePoint(A, B);
            Vector3D v = Vector3D.CrossProduct(A.Normal, B.Normal);
            return new LineVector(v, aPoint);
        }


        public static Point3D SolveForOnePoint(Plane A, Plane B)
        {
            Vector othertwopoints = new Vector();
            Point3D p = new Point3D();
            Vector[] ret = Return2SVectorGivenOmissionOfElementFromNormal("X", A, B);
            if (IsDeterminantZero(ret) == false)
            {
                othertwopoints = SolveSquareMatrixGiven3Vectors(ret);
                p.X = 0;
                p.Y = othertwopoints.X;
                p.Z = othertwopoints.Y;
                return p;
            }

            ret = Return2SVectorGivenOmissionOfElementFromNormal("Y", A, B);
            if (IsDeterminantZero(ret) == false)
            {
                othertwopoints = SolveSquareMatrixGiven3Vectors(ret);
                p.Y = 0;
                p.X = othertwopoints.X;
                p.Z = othertwopoints.Y;

                return p;
            }
            ret = Return2SVectorGivenOmissionOfElementFromNormal("Z", A, B);
            othertwopoints = SolveSquareMatrixGiven3Vectors(ret);
            p.X = othertwopoints.X;
            p.Y = othertwopoints.Y;
            p.Z = 0;
            return p;


        }
        public static Vector SolveSquareMatrixGiven3Vectors(Vector[] cols)
        {
            double a1 = 0;
            double a2 = 0;
            double b1 = 0;
            double b2 = 0;
            double d1 = 0;
            double d2 = 0;
            a1 = cols[0].X;
            a2 = cols[0].Y;
            b1 = cols[1].X;
            b2 = cols[1].Y;
            d1 = cols[2].X;
            d2 = cols[2].Y;
            double denom = b2 * a1 - b1 * a2;
            double x1 = (b2 * d1 - b1 * d2) / denom;
            double x2 = (d2 * a1 - d1 * a2) / denom;
            return new Vector(x1, x2);

        }
        public static bool IsDeterminantZero(Vector[] cols)
        {
            //// a  b
            ////  c  d
            /// above is the matrix
            /// 
            double a = cols[0].X;
            double c = cols[0].Y;

            double b = cols[1].X;
            double d = cols[1].Y;

            double determinant = a * d - b * c;
            return (determinant == 0);
        }
        public static Vector[] Return2SVectorGivenOmissionOfElementFromNormal(string WHatToOmit, Plane A, Plane B)
        {
            Vector col1 = new Vector();
            Vector col2 = new Vector();
            Vector col3 = new Vector();
            switch (WHatToOmit)
            {
                case "X":
                    col1.X = A.Normal.Y;
                    col1.Y = B.Normal.Y;
                    col2.X = A.Normal.Z;
                    col2.Y = B.Normal.Z;
                    break;
                case "Y":
                    col1.X = A.Normal.X;
                    col1.Y = B.Normal.X;
                    col2.X = A.Normal.Y;
                    col2.Y = B.Normal.Y;
                    break;
                case "Z":
                    col1.X = A.Normal.X;
                    col1.Y = B.Normal.X;
                    col2.X = A.Normal.Y;
                    col2.Y = B.Normal.Y;
                    break;
                default:
                    break;
            }
            col3.X = A.Equation.Constantparameter;
            col3.Y = B.Equation.Constantparameter;
            Vector[] ret = new Vector[] { col1, col2, col3 };
            return ret;
            //Matrix ret = new Matrix();
            //ret.
        }

        public Plane(Vector3D aNormal, Point3D aPoint)
        {
            Normal = aNormal;
            PointOnPlane = aPoint;
            Equation = new PlaneEquation(this);
        }

    }
    public class LineVector
    {
        Vector3D prVectorDirection;
        Point3D prPointOnPlane;
        public Vector3D VectorDirection
        {
            get
            {
                return prVectorDirection;
            }
            set
            {
                prVectorDirection = value;
            }
        }
        public Point3D PointOnPlan
        {
            get
            {
                return prPointOnPlane;
            }
            set
            {
                prPointOnPlane = value;
               
            }
        }
        public LineVector(Vector3D dir,Point3D p)
        {
            VectorDirection = dir;
            PointOnPlan = p;

        }
        public LineVector()
        {

        }
    }
    public class PlaneEquation
    {

        Plane plane;
        double  prXparameter;
        double prYparameter;
        double prZparameter;
        double prConstantparameter;
        // asumption is that the quation is of the form a1*X+a2*y+a3*z=C;
        public double Xparameter
        {
            get
            {
                return prXparameter;
            }
            set
            {
                prXparameter = value;
            }
        }
        public double Yparameter
        {
            get
            {
                return prYparameter;
            }
            set
            {
                prYparameter = value;
            }
        }
        public double Zparameter
        {
            get
            {
                return prZparameter;
            }
            set
            {
                prZparameter = value;
            }
        }
        public double Constantparameter
        {
            get
            {
                return prConstantparameter;
            }
            set
            {
                prConstantparameter = value;
            }
        }
        
        public PlaneEquation()
        {

        }
        
        public PlaneEquation(Plane aplane)
        {
            plane = aplane;
            Vector3D normal=plane.Normal;
            Vector3D poinonPlane=(Vector3D)plane.PointOnPlane;
            Xparameter = normal.X;
            Yparameter = normal.Y;
            Zparameter = normal.Z;
            Constantparameter = Vector3D.DotProduct(normal, poinonPlane);
        }
        public bool solveForRemainingVariable (EquationParameter EP,out double solvedForVAlue) 
        {
            solvedForVAlue = 0;
            bool result;
            double ret = 0;
            double denominator = 0;
            double numerator = 0;
            if (EP.X == null)
            {
                numerator = prConstantparameter - Yparameter * (double)EP.Y - Zparameter * (double)EP.Z;
                denominator = Xparameter;
            }
            else if (EP.Y == null)
            {
                numerator = prConstantparameter - Zparameter * (double)EP.Z - Xparameter * (double)EP.X;
                denominator = Yparameter;
            }
            else if (EP.Z == null)
            {
                numerator = prConstantparameter - Yparameter * (double)EP.Y - Xparameter * (double)EP.X;
                denominator = Zparameter;
            }
            if (denominator == 0)
            {
                result = false;
            }
            else
            {
                result = true;
                solvedForVAlue = numerator / denominator;
            }
            return result;
        }
    }
}
public class EquationParameter
{
    private object x;
    private object y;
    private object z;
    public object X
    {
        get
        {
            return x;
        }
        set
        {
            x = value;
        }
    }
    public object Z
    {
        get
        {
            return z;
        }
        set
        {
            z = value;
        }
    }
    public object Y
    {
        get
        {
            return y;
        }
        set
        {
            y = value;
        }
    }
    public void RotatePrameters()
    {
        // x==> Y & Y ==> Z & Z ==> X
        object temp = Z;
        Z = X;
        X = Y;
        Y = temp;
    }
}
