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
    public class VectorProjections
    {
        public static Vector3D calculateProjectionOfVectorIntoPlane(Vector3D A,Vector3D planeNormal,out Vector3D distanceVector)
        {
            Vector3D[] basisOfPlane = ReturnPerpendicularOrthognalBasisPerNormal(planeNormal);
            Vector3D projection = basisOfPlane[0] * (Vector3D.DotProduct(A, basisOfPlane[0]) / Vector3D.DotProduct(basisOfPlane[0], basisOfPlane[0])) + basisOfPlane[1] * (Vector3D.DotProduct(A, basisOfPlane[1]) / Vector3D.DotProduct(basisOfPlane[1], basisOfPlane[1]));
            distanceVector=A-projection;
            return projection;
        }

        public static Vector3D[] ReturnPerpendicularOrthognalBasisPerNormal(Vector3D NormalVector)
        {
            Vector3D orthognalFirst = CalculateGroundAxixGivenTubePerpendicularAxis(NormalVector);
            orthognalFirst = orthognalFirst / orthognalFirst.Length;
            Vector3D orthognalSecond = Vector3D.CrossProduct(NormalVector / NormalVector.Length, orthognalFirst);


            return new Vector3D[2] { orthognalFirst, orthognalSecond };
        }
        public static Vector3D CalculateGroundAxixGivenTubePerpendicularAxis(Vector3D NormalVector)
        {

            double RetVectorX;
            double RetVectorY;
            double RetVectorZ;
            double PerpenX = NormalVector.X;
            double PerpenZ = NormalVector.Z;
            double PerpenY = NormalVector.Y;
            if (!(PerpenZ == 0))
            {

                RetVectorX = 1;
                RetVectorY = 0;//  we are on the ground plane

                double dotProductMultiplication_knowPart = PerpenX * RetVectorX + RetVectorY * RetVectorY;
                RetVectorZ = -dotProductMultiplication_knowPart / PerpenZ;
            }
            else if (!(PerpenY == 0))
            {
                RetVectorX = 1;
                RetVectorZ = 0;//  we are on the ground plane

                double dotProductMultiplication_knowPart = PerpenX * RetVectorX + RetVectorZ * RetVectorZ;
                RetVectorY = -dotProductMultiplication_knowPart / PerpenY;
            }
            else
            {
                RetVectorY = 1;
                RetVectorZ = 0;//  we are on the ground plane

                double dotProductMultiplication_knowPart = PerpenY * RetVectorY + RetVectorZ * RetVectorZ;
                RetVectorX = -dotProductMultiplication_knowPart / PerpenX;
            }
            Vector3D ret = new Vector3D(RetVectorX, RetVectorY, RetVectorZ);
            return ret;

        }

        public static bool findIntersectionOfVectorOnPalneStatringAtPointP(Point3D p, Vector3D v, Vector3D Normal_orig, out Point3D final)
        {
            bool result;
            Vector3D Normal=Normal_orig;
            Normal.Normalize();
            /// s is the sought length of our vector v; the right length will make it just as long to touch the plane
            // (p -sV) dot n=0 ==> s=(p dot n)/(v dot n)
            double dotProductvdotn = Vector3D.DotProduct(v, Normal);
            double dotProductpdotn = Vector3D.DotProduct((Vector3D)p, Normal);
            if (dotProductvdotn == 0)
            {
                result = false;
                final = new Point3D();
            }
            else
            {
                result = true;
                double s = dotProductpdotn / dotProductvdotn;
                final = p - s * v;
            }

            return result;

        }
    }
}
