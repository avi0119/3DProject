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
    public class Cube_old : Primitive
    {
        Model3DGroup prGroup;
        public Model3DGroup group
        {
            get;
            set;
        }
        PointCollection textCoord_alternative = new PointCollection();
        Point3D[] positions = new Point3D[] {
        new Point3D(0, 1, 1),
        new Point3D(0, 0, 1),
        new Point3D(1, 1, 1),
        new Point3D(1, 0, 1),
        new Point3D(1, 1, 0),
        new Point3D(0, 0, 0),
        new Point3D(0, 1, 1),
        new Point3D(0, 0, 1),
        new Point3D(0, 1, 0),
        new Point3D(0, 1, 0),
        new Point3D(1, 1, 1),
        new Point3D(0, 1, 1),
        new Point3D(1, 1, 1),
        new Point3D(1, 0, 1),
        new Point3D(1, 1, 0),
        new Point3D(1, 0, 0),
        new Point3D(1, 0, 1),
        new Point3D(0, 0, 1),
        new Point3D(1, 0, 0),
        new Point3D(0, 0, 0),
        new Point3D(1, 1, 0),
        new Point3D(1, 0, 0),
        new Point3D(0, 1, 0),
        new Point3D(0, 0, 0)
        };
        int[] vertices = new int[] {
        0,
        1,
        2,

        1,
        3,
        2,

        4,
        5,
        6,

        5,
        7,
        6,

        8,
        9,
        10,

        9,
        11,
        10,

        12,
        13,
        14,

        13,
        15,
        14,

        16,
        17,
        18,

        17,
        19,
        18,

        20,
        21,
        22,

        21,
        23,
        22
        };
        private void AddOneFramePer6Positins()
        {
            textCoord_alternative.Add(new Point(0,0));
            textCoord_alternative.Add(new Point(1,0));
            textCoord_alternative.Add(new Point(1,1));
            textCoord_alternative.Add(new Point(1,0));
            textCoord_alternative.Add(new Point(0,1));
            textCoord_alternative.Add(new Point(1,1));

            textCoord_alternative.Add(new Point(0, 0));
            textCoord_alternative.Add(new Point(1, 0));
            textCoord_alternative.Add(new Point(1, 1));
            textCoord_alternative.Add(new Point(1, 0));
            textCoord_alternative.Add(new Point(0, 1));
            textCoord_alternative.Add(new Point(1, 1));

            textCoord_alternative.Add(new Point(0, 0));
            textCoord_alternative.Add(new Point(1, 0));
            textCoord_alternative.Add(new Point(1, 1));
            textCoord_alternative.Add(new Point(1, 0));
            textCoord_alternative.Add(new Point(0, 1));
            textCoord_alternative.Add(new Point(1, 1));

            textCoord_alternative.Add(new Point(0, 0));
            textCoord_alternative.Add(new Point(1, 0));
            textCoord_alternative.Add(new Point(1, 1));
            textCoord_alternative.Add(new Point(1, 0));
            textCoord_alternative.Add(new Point(0, 1));
            textCoord_alternative.Add(new Point(1, 1));

            textCoord_alternative.Add(new Point(0, 0));
            textCoord_alternative.Add(new Point(1, 0));
            textCoord_alternative.Add(new Point(1, 1));
            textCoord_alternative.Add(new Point(1, 0));
            textCoord_alternative.Add(new Point(0, 1));
            textCoord_alternative.Add(new Point(1, 1));

            textCoord_alternative.Add(new Point(0, 0));
            textCoord_alternative.Add(new Point(1, 0));
            textCoord_alternative.Add(new Point(1, 1));
            textCoord_alternative.Add(new Point(1, 0));
            textCoord_alternative.Add(new Point(0, 1));
            textCoord_alternative.Add(new Point(1, 1));

            textCoord_alternative.Add(new Point(0, 0));
            textCoord_alternative.Add(new Point(1, 0));
            textCoord_alternative.Add(new Point(1, 1));
            textCoord_alternative.Add(new Point(1, 0));
            textCoord_alternative.Add(new Point(0, 1));
            textCoord_alternative.Add(new Point(1, 1));


        }
        public Cube_old()
        {
            AddOneFramePer6Positins();
            this.Mesh.Positions = new Point3DCollection(positions);
            this.Mesh.TriangleIndices = new Int32Collection(vertices);
            this.Mesh.TextureCoordinates = textCoord_alternative;
            object ib = TubeGeometry.CreateTransparemtFramedBrush();

            //MeshGeometry3D mesh = this.Mesh;
            //MaterialGroup MatGR = new MaterialGroup();
            //MaterialGroup CircelMaterialGroup = new MaterialGroup();
            //MaterialGroup AxisPoleMaterialGroup = new MaterialGroup();
            //MaterialGroup BAckMatGR = new MaterialGroup();
            DiffuseMaterial m2;
            if (ib is ImageBrush)
            {
                //MatGR.Children.Add(new DiffuseMaterial(ib as ImageBrush));
                m2 = new DiffuseMaterial(ib as ImageBrush);
            }
            else
            {
                //MatGR.Children.Add(new DiffuseMaterial(ib as DrawingBrush));
                m2 = new DiffuseMaterial(ib as DrawingBrush);
            }


            
  







            Model3DGroup aGroup = new Model3DGroup();
            var material = new DiffuseMaterial(System.Windows.Media.Brushes.Green);
            material = m2;
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
    }

    public abstract class Primitive
    {
        public Primitive()
        {
            this.Mesh = new MeshGeometry3D();
        }
        public MeshGeometry3D Mesh { get; protected set; }
    }
    public class Cube : Primitive
    {
        Model3DGroup prGroup;
        public Model3DGroup group
        {
            get;
            set;
        }
        PointCollection textCoord_alternative = new PointCollection();
        Point3D[] positions = new Point3D[] {

                /// top
        
        new Point3D(0, 1, 1),
        new Point3D(1, 1, 1),
        new Point3D(1, 1, 0),
        new Point3D(0, 1, 0),


        //east
        new Point3D(1, 0, 1),
        new Point3D(1, 0, 0),
        new Point3D(1, 1, 0),
        new Point3D(1, 1, 1)
        ,


                //west
        new Point3D(0, 0, 1),
        new Point3D(0, 0, 0),
        new Point3D(0, 1, 0),
        new Point3D(0, 1, 1),

                
        //north
        new Point3D(0, 0, 0),
        new Point3D(1, 0, 0),
        new Point3D(1, 1, 0),
        new Point3D(0, 1, 0),


        //south

        new Point3D(0, 0, 1),
        new Point3D(1, 0, 1),
        new Point3D(1, 1, 1),
        new Point3D(0, 1, 1)
        ,

        // bottom
        new Point3D(0, 0, 0),
        new Point3D(1, 0, 0),
        new Point3D(1, 0, 1),
        new Point3D(0, 0, 1)
        









        };
        int[] vertices = new int[] {
            // top
        0,
        1,
        2,

        0,
        2,
        3
       
        ,

        // east


        4,
        5,
        6,

        4,
        6,
        7
        ,

                ///west
        
        8,
        9,
        10,

        8,
        10,
        11
        ,
        // north

        12,
        13,
        14,

        12,
        14,
        15
        ,

        // south

        16,
        17,
        18,

        16,
        18,
        19
        ,

         
        // bottom
        20,
        21,
        22,

        20,
        22,
        23

        };
        private void AddOneFramePer6Positins()
        {
            for (int i = 0; i < 6; i++)
            {
                textCoord_alternative.Add(new Point(0, 0));
                textCoord_alternative.Add(new Point(1, 0));
                textCoord_alternative.Add(new Point(1, 1));
                textCoord_alternative.Add(new Point(0, 1));
                //textCoord_alternative.Add(new Point(1, 1));
                //textCoord_alternative.Add(new Point(0, 1));
            }

  


        }
        public Cube(bool increase =false)
        {
            if (increase==true)
            {
                var k = positions.Select(i => { Vector3D vv = (Vector3D)i * 100; return (Point3D)vv; });
                positions = k.ToArray();
            }
            AddOneFramePer6Positins();
            this.Mesh.Positions = new Point3DCollection(positions);
            this.Mesh.TriangleIndices = new Int32Collection(vertices);

            this.Mesh.TextureCoordinates = textCoord_alternative;

            object ib = TubeGeometry.CreateTransparemtFramedBrush();

            //MeshGeometry3D mesh = this.Mesh;
            //MaterialGroup MatGR = new MaterialGroup();
            //MaterialGroup CircelMaterialGroup = new MaterialGroup();
            //MaterialGroup AxisPoleMaterialGroup = new MaterialGroup();
            //MaterialGroup BAckMatGR = new MaterialGroup();
            DiffuseMaterial m2;
            if (ib is ImageBrush)
            {
                //MatGR.Children.Add(new DiffuseMaterial(ib as ImageBrush));
                m2 = new DiffuseMaterial(ib as ImageBrush);
            }
            else
            {
                //MatGR.Children.Add(new DiffuseMaterial(ib as DrawingBrush));
                m2 = new DiffuseMaterial(ib as DrawingBrush);
            }











            Model3DGroup aGroup = new Model3DGroup();
            var material = new DiffuseMaterial(System.Windows.Media.Brushes.Green);

            material = m2;

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
    }
}
