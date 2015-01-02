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
    public class DrawBar
    {
        public Model3DGroup group
        {
            get;
            set;
        }
        MeshGeometry3D prMesh;
        public MeshGeometry3D Mesh
        {
            get
            {
                return prMesh;
            }
            set
            {
                prMesh = value;
            }
        }
        public DrawBar()
        {

        }
        public DrawBar(double step, Point3D A, Point3D B,double barDiameter, Brush color)
        {
            Vector3D NormalVector = A - B;
            NormalVector.Normalize();
            Vector3D[] basis = ReturnPerpendicularOrthognalBasisPerNormal(NormalVector);
            Point3D[] arr = ReturnArrayOfPointAlongTheCircumferenceGivenOrthignalBasis(basis,step,barDiameter/2);
            int arrCount = arr.Length;

            MeshGeometry3D mesh = new MeshGeometry3D();
            Point3DCollection P3DcollectionPositionsOneSide = returnPointsOnCircomferenceGivenTranslationVector(arr,(Vector3D)A);
            Point3DCollection P3DcollectionPositionsSecondSide = returnPointsOnCircomferenceGivenTranslationVector(arr, (Vector3D)B);
            var res=P3DcollectionPositionsOneSide.Concat(P3DcollectionPositionsSecondSide);

           
            Point3D[] arrayOf3DPoints = res.ToArray();
            Point3DCollection finalpPoint3DCollection = new Point3DCollection(arrayOf3DPoints);

            mesh.Positions = finalpPoint3DCollection;
            Int32Collection Triangles_int = moveCrissCrossFromOnePointToAnotherAndStoreMovemenstIntoIntCollection(arrCount);

            mesh.TriangleIndices = Triangles_int;


            //mesh.Normals = normals;

            //mesh.TextureCoordinates = textCoord_alternative;
            Mesh = mesh;


            Model3DGroup aGroup = new Model3DGroup();
            var material = new DiffuseMaterial(color);

            //material = m2;

            //MeshGeometry3D bodyMesh = CreateMeshArrowBody();
            //MeshGeometry3D headMesh = CreateMeshArrowHead();
            GeometryModel3D a;
            a = new GeometryModel3D(this.Mesh, material);
            a.BackMaterial = material;
            aGroup.Children.Add(a);


            //GeometryModel3D b;
            //material = new DiffuseMaterial(System.Windows.Media.Brushes.Yellow);
            //b = new GeometryModel3D(headMesh, material);
            //b.BackMaterial = material;
            //aGroup.Children.Add(b);
            group = aGroup;

        }
        public Int32Collection  moveCrissCrossFromOnePointToAnotherAndStoreMovemenstIntoIntCollection(int countOfPointsInEachEdge)
        {

            Int32Collection ret = new Int32Collection();
            for (int counter = 0; counter < countOfPointsInEachEdge; counter++)
            {
                if (counter + 1 < countOfPointsInEachEdge)
                {
                    // normal criss crossing 
                    ret.Add(counter);
                    ret.Add(counter + countOfPointsInEachEdge);
                    ret.Add(counter+1);

                    ret.Add(counter + countOfPointsInEachEdge);
                    ret.Add(counter +1+ countOfPointsInEachEdge);
                    ret.Add(counter + 1);

                }
                else
                {
                    // last vertical  line need to connect to first line
                    ret.Add(counter);
                    ret.Add(counter + countOfPointsInEachEdge);
                    ret.Add(0);

                    ret.Add(counter + countOfPointsInEachEdge);
                    ret.Add( countOfPointsInEachEdge);
                    ret.Add(0);
                }
            }
            return ret;
        }
        public Point3DCollection returnPointsOnCircomferenceGivenTranslationVector(Point3D[] arr,Vector3D translationVect)
        {
            var result = arr.Select(element => translationVect+element);
            return new Point3DCollection(result);
        }
        public Point3D[] ReturnArrayOfPointAlongTheCircumferenceGivenOrthignalBasis(Vector3D[] basis,double step,double radius)
        {
            List<Point3D> listPoints = new List<Point3D>();
            double sin = 0; double cos = 0;
            double pi = Math.PI;
            Vector3D currVectorRepresentingPoint = new Vector3D();
            Point3D currPoint = new Point3D();
            Vector3D v1= basis[0];
            Vector3D v2= basis[1];
            for (double angle = 0; angle < 360; angle = angle + step)
            {
                sin = Math.Sin(angle / 180 * pi);
                cos = Math.Cos(angle / 180 * pi);
                currVectorRepresentingPoint = (v1 * sin  + v2 * cos )* radius;
                currPoint = (Point3D)currVectorRepresentingPoint;
                listPoints.Add(currPoint);

            }
            return  listPoints.ToArray();
        }
        public Vector3D[] ReturnPerpendicularOrthognalBasisPerNormal(Vector3D NormalVector)
        {
            Vector3D orthognalFirst = CalculateGroundAxixGivenTubePerpendicularAxis(NormalVector);
            orthognalFirst = orthognalFirst / orthognalFirst.Length;
            Vector3D orthognalSecond = Vector3D.CrossProduct(NormalVector / NormalVector.Length, orthognalFirst);


            return new Vector3D[2] { orthognalFirst, orthognalSecond };
        }
        public Vector3D CalculateGroundAxixGivenTubePerpendicularAxis(Vector3D NormalVector)
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
    }
}
