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
    class Arrow3DShape
    {
        Model3DGroup prGroup;
        public Model3DGroup group
        {
            get;
            set;
        }
        public Arrow3DShape()
        {
            Model3DGroup aGroup = new Model3DGroup();
            var material = new DiffuseMaterial(System.Windows.Media.Brushes.Green);
            MeshGeometry3D bodyMesh = CreateMeshArrowBody();
            MeshGeometry3D headMesh = CreateMeshArrowHead();
            GeometryModel3D a;
            a = new GeometryModel3D(bodyMesh, material);
            a.BackMaterial = material;
            aGroup.Children.Add(a);


            GeometryModel3D b;
             material = new DiffuseMaterial(System.Windows.Media.Brushes.Yellow);
            b = new GeometryModel3D(headMesh, material);
            b.BackMaterial = material;
            aGroup.Children.Add(b);
            group = aGroup;

        }

        private MeshGeometry3D CreateMeshArrowBody()
        {

            Vector3DCollection normals = new Vector3DCollection();
            Point3DCollection P3DcollectionPositions = new Point3DCollection();
            Int32Collection Triangles_int = new Int32Collection();
            PointCollection textCoord = new PointCollection();
            PointCollection textCoord_alternative = new PointCollection();
            MeshGeometry3D triangleMesh = new MeshGeometry3D();
            double bodyheight = 6;
            P3DcollectionPositions.Add(new Point3D(1,0,1)); // 0
            P3DcollectionPositions.Add(new Point3D(1, 0, -1)); // 1
            P3DcollectionPositions.Add(new Point3D(-1, 0, -1)); // 2
            P3DcollectionPositions.Add(new Point3D(-1, 0, 1)); // 3 

            P3DcollectionPositions.Add(new Point3D(1, bodyheight, 1)); // 4
            P3DcollectionPositions.Add(new Point3D(1, bodyheight, -1)); // 5
            P3DcollectionPositions.Add(new Point3D(-1, bodyheight, -1)); // 6
            P3DcollectionPositions.Add(new Point3D(-1, bodyheight, 1)); // 7 


            // arrow bottom

            Triangles_int.Add(0);
            Triangles_int.Add(1);
            Triangles_int.Add(2);

            Triangles_int.Add(0);
            Triangles_int.Add(2);
            Triangles_int.Add(3);

            // arrow first side

            Triangles_int.Add(0);
            Triangles_int.Add(4);
            Triangles_int.Add(1);

            Triangles_int.Add(1);
            Triangles_int.Add(4);
            Triangles_int.Add(5);

            // arrow second side
            Triangles_int.Add(0);
            Triangles_int.Add(4);
            Triangles_int.Add(3);

            Triangles_int.Add(4);
            Triangles_int.Add(7);
            Triangles_int.Add(3);


            //// arrow third side

            Triangles_int.Add(3);
            Triangles_int.Add(2);
            Triangles_int.Add(6);

            Triangles_int.Add(3);
            Triangles_int.Add(7);
            Triangles_int.Add(2);

            //// arrow fourth side

            Triangles_int.Add(2);
            Triangles_int.Add(1);
            Triangles_int.Add(5);

            Triangles_int.Add(5);
            Triangles_int.Add(6);
            Triangles_int.Add(2);

            // arrow top

            Triangles_int.Add(4);
            Triangles_int.Add(5);
            Triangles_int.Add(6);

            Triangles_int.Add(4);
            Triangles_int.Add(6);
            Triangles_int.Add(7);



            triangleMesh.Positions = P3DcollectionPositions;
            //triangleMesh.Normals = normals;
            triangleMesh.TriangleIndices = Triangles_int;
            return triangleMesh;
        }
        private MeshGeometry3D CreateMeshArrowHead()
        {

            Vector3DCollection normals = new Vector3DCollection();
            Point3DCollection P3DcollectionPositions = new Point3DCollection();
            Int32Collection Triangles_int = new Int32Collection();
            PointCollection textCoord = new PointCollection();
            PointCollection textCoord_alternative = new PointCollection();
            MeshGeometry3D triangleMesh = new MeshGeometry3D();
            double headhight = 3;
            double headbottomsurfacehaldside = 1.5;
            double bodyheight = 6;
            P3DcollectionPositions.Add(new Point3D(headbottomsurfacehaldside, bodyheight, headbottomsurfacehaldside)); // 0
            P3DcollectionPositions.Add(new Point3D(headbottomsurfacehaldside, bodyheight, -headbottomsurfacehaldside)); // 1
            P3DcollectionPositions.Add(new Point3D(-headbottomsurfacehaldside, bodyheight, -headbottomsurfacehaldside)); // 2
            P3DcollectionPositions.Add(new Point3D(-headbottomsurfacehaldside, bodyheight, headbottomsurfacehaldside)); // 3 
            P3DcollectionPositions.Add(new Point3D(0, bodyheight+headhight, 0)); // 4



            // arrow bottom

            Triangles_int.Add(0);
            Triangles_int.Add(1);
            Triangles_int.Add(2);

            Triangles_int.Add(0);
            Triangles_int.Add(2);
            Triangles_int.Add(3);

            // arrow first side

            Triangles_int.Add(0);
            Triangles_int.Add(1);
            Triangles_int.Add(4);

            // arrow second side

            Triangles_int.Add(1);
            Triangles_int.Add(2);
            Triangles_int.Add(4);

            // arrow third side

            Triangles_int.Add(2);
            Triangles_int.Add(3);
            Triangles_int.Add(4);



            // arrow fourth side

            Triangles_int.Add(3);
            Triangles_int.Add(0);
            Triangles_int.Add(4);



            triangleMesh.Positions = P3DcollectionPositions;
            //triangleMesh.Normals = normals;
            triangleMesh.TriangleIndices = Triangles_int;
            return triangleMesh;
        }
    }
}
