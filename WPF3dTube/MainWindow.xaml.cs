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
using System.Drawing;
using System.Resources;

namespace WPF3dTube
{

    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public AxisAngleRotation3D rotationY = new AxisAngleRotation3D();
        public MatrixTransform3D myMatrixTransform3D;
        private GeometryModel3D mGeometry;
        static TimeSpan lastRenderTime = new TimeSpan();
        static private double selfrotatingangle=0;
        static private double groundrotatingangle = 0;
        private int SpinMode = 0;

        public MainWindow()
        {
            InitializeComponent();
            //UpdateLayout();
            Build3DTube();
            CompositionTarget.Rendering += CompositionTarget_Rendering;
            
            
        }
        void CompositionTarget_Rendering(object sender, EventArgs e)
        {
            //return;
            // Ensure we only do this once per frame
            Vector3D TestVector = new Vector3D(1, 1, 1);
            Vector3D RestlVectoTest;
            Vector3D RestlVectoTest2;
            Vector3D RestlVectoTest3;
            if (lastRenderTime == ((RenderingEventArgs)e).RenderingTime)
                return;

            lastRenderTime = ((RenderingEventArgs)e).RenderingTime;
            //rotationY.Angle = rotationY.Angle + 1;
            //if (SpinMode <= 2)
            //{
            //    selfrotatingangle = selfrotatingangle + 2;
            //}
            //else
            //{
            //    groundrotatingangle = groundrotatingangle + 2;
            //}

                selfrotatingangle = selfrotatingangle + 1.6;

                groundrotatingangle = groundrotatingangle + .8;
 
            RecalculateTheAxis(groundrotatingangle);
            Matrix3D matrix1 = returnMatrixPertransformation(selfrotatingangle, TubeGeometry.TheAxis);//selfrotatingangle was removed
            Matrix3D matrix2 = returnMatrixPertransformation(groundrotatingangle, TubeGeometry.GroundAxis);//groundrotatingangle  was removed
            RestlVectoTest = Vector3D.Multiply(TestVector, matrix1);
            RestlVectoTest2 = Vector3D.Multiply(RestlVectoTest, matrix2); ;
            //Matrix3D retmtarix3D = Matrix3D.Multiply(matrix1, matrix2);
            //matrix2.Append(matrix1);
            matrix1.Append(matrix2);
            RestlVectoTest3 = Vector3D.Multiply(TestVector, matrix1); ;
            //mGeometry.Transform.Value = matrix2;
            //myMatrixTransform3D = matrix2;
            ((Transform3DGroup)(mGeometry.Transform)).Children[1] = new MatrixTransform3D(matrix1);

            SpinMode = SpinMode + 1;
            if (SpinMode == 4) SpinMode = 0;


        }
        void CompositionTarget_Rendering_old(object sender, EventArgs e)
        {
            //return;
            // Ensure we only do this once per frame
            Vector3D TestVector = new Vector3D(1, 1, 1);
            Vector3D RestlVectoTest;
            Vector3D RestlVectoTest2;
            Vector3D RestlVectoTest3;
            if (lastRenderTime == ((RenderingEventArgs)e).RenderingTime)
                return;

            lastRenderTime = ((RenderingEventArgs)e).RenderingTime;
            //rotationY.Angle = rotationY.Angle + 1;
            if (SpinMode <= 2)
            {
                selfrotatingangle = selfrotatingangle + .1;
            }
            else
            {
                groundrotatingangle = groundrotatingangle + .5;
            }
            RecalculateTheAxis(groundrotatingangle);
            Matrix3D matrix1 = returnMatrixPertransformation(selfrotatingangle, TubeGeometry.TheAxis);//selfrotatingangle was removed
            Matrix3D matrix2 = returnMatrixPertransformation(groundrotatingangle, TubeGeometry.GroundAxis);//groundrotatingangle  was removed
            RestlVectoTest = Vector3D.Multiply(TestVector, matrix1);
            RestlVectoTest2 = Vector3D.Multiply(RestlVectoTest, matrix2); ;
            //Matrix3D retmtarix3D = Matrix3D.Multiply(matrix1, matrix2);
            //matrix2.Append(matrix1);
            matrix1.Append(matrix2);
            RestlVectoTest3 = Vector3D.Multiply(TestVector, matrix1); ;
            //mGeometry.Transform.Value = matrix2;
            //myMatrixTransform3D = matrix2;
            ((Transform3DGroup)(mGeometry.Transform)).Children[1] = new MatrixTransform3D(matrix1);

            SpinMode = SpinMode + 1;
            if (SpinMode == 4) SpinMode = 0;


        }
        private void RecalculateTheAxis(double groundrotatingangle)
        { 
            /// the tube rotated around GroundAxis by  degrees
            /// therefore it is necessary to rotate OriginalAxis with respect to GroundAxis By groundrotatingangle
            ///  that is done by setting TheAxis to the OriginalAxis rotated by groundrotatingangle
            
        }
        private Matrix3D returnMatrixPertransformation(double angle,Vector3D TheAxis )
        {
            Transform3DGroup Tr3D = new Transform3DGroup();
            RotateTransform3D myRotateTransform3D = new RotateTransform3D();
            AxisAngleRotation3D myAxisAngleRotation3d = new AxisAngleRotation3D();
            myAxisAngleRotation3d.Axis = TheAxis;
            myAxisAngleRotation3d.Angle = angle;
            myRotateTransform3D.Rotation = myAxisAngleRotation3d;

            Tr3D.Children.Add(myRotateTransform3D);
            //mGeometry.Transform = Tr3D;
            return Tr3D.Value;
        }
        void CompositionTarget_Rendering_old2(object sender, EventArgs e)
        {
            //return;
            // Ensure we only do this once per frame
            if (lastRenderTime == ((RenderingEventArgs)e).RenderingTime)
                return;

            lastRenderTime = ((RenderingEventArgs)e).RenderingTime;
            rotationY.Angle = rotationY.Angle + 1;

            Matrix3D Current = mGeometry.Transform.Value;
            groundrotatingangle = groundrotatingangle + 1;
            Current.Rotate(new Quaternion(TubeGeometry.GroundAxis, groundrotatingangle));

        }
        public void Build3DTube()
        {

             

            ImageBrush ib = new ImageBrush();

            //ib.ImageSource = new BitmapImage(new Uri("IMG_0804.jpg", UriKind.Relative));
            //ib.ImageSource = new BitmapImage(new Uri(@"Resources\IMG_0804.JPG", UriKind.Relative));
            //ib.ImageSource = new BitmapImage(new Uri("/Resources/IMG_0804.JPG", UriKind.Relative));
            Bitmap ResourceBitmap = (Bitmap)WPF3dTube.Properties.Resources.IMG_0804;

            BitmapSource bs = System.Windows.Interop.Imaging.CreateBitmapSourceFromHBitmap(
                                              ResourceBitmap.GetHbitmap(),
                                              IntPtr.Zero,
                                              System.Windows.Int32Rect.Empty,
                                              BitmapSizeOptions.FromWidthAndHeight(ResourceBitmap.Width, ResourceBitmap.Height));
            ImageBrush ib2 = new ImageBrush(bs);


            ib.ImageSource = new BitmapImage(new Uri(@"IMG_0804.JPG", UriKind.Relative));
            ib = ib2;
            //ib.ImageSource = rm;
            
            //Bitmap myImage = (Bitmap)rm.GetObject("myImage");
            //Bitmap bmp = new Bitmap(System.Reflection.Assembly.GetEntryAssembly().GetManifestResourceStream("WPF3dTube.Resources.IMG_0804.JPG"));
            //Bitmap myImage = (Bitmap)WPF3dTube.Properties.Resources.ResourceManager.GetObject("desert");
            //string[] all = System.Reflection.Assembly.GetEntryAssembly().GetManifestResourceNames();

            //foreach (string one in all)
            //{
            //    MessageBox.Show(one);
            //}
            Vector3DCollection normals = new Vector3DCollection();
            Point3DCollection P3DcollectionPositions = new Point3DCollection();
            Int32Collection Triangles_int = new Int32Collection();
            PointCollection textCoord = new PointCollection();
            PointCollection textCoord_alternative = new PointCollection();
            Vector3D AxisVector;
            AxisVector = new Vector3D(1, 1, 0);
            TubeGeometry.GroundAxis = TubeGeometry.FindGroundAxixGivenTubePerpendicularAxis(AxisVector);
            double TubeThickness = 3;
            double InternalDiameter = 8;
            int angleStep_horizotal = 2;
            int angleStep_vertical = 2;
            int NumberOfDefreesToComplete_horizontal = 360;
            int NumberOfDefreesToComplete_vertical = 360;
            TubeGeometry.ReturnPositionsTriangleIndicesNormals(NumberOfDefreesToComplete_horizontal, NumberOfDefreesToComplete_vertical, angleStep_horizotal,angleStep_vertical, TubeThickness, InternalDiameter, AxisVector, ref  P3DcollectionPositions, ref  Triangles_int, ref  normals, ref textCoord,ref  textCoord_alternative);
            MeshGeometry3D mesh = new MeshGeometry3D();

            mesh.Positions = P3DcollectionPositions;
            mesh.Normals = normals;
            mesh.TriangleIndices = Triangles_int;
            //mesh.TextureCoordinates = textCoord;
            mesh.TextureCoordinates = textCoord_alternative;
            MaterialGroup MatGR = new MaterialGroup();
            MaterialGroup BAckMatGR = new MaterialGroup();
            //MatGR.Children.Add(new DiffuseMaterial(TubeGeometry.ReturnLinearBrush()));
            //MatGR.Children.Add(new DiffuseMaterial(TubeGeometry.CreateTransparemtFramedBrush())); /// commented wire-like
            MatGR.Children.Add(new DiffuseMaterial(ib));
            //MatGR.Children.Add(new DiffuseMaterial(TubeGeometry.CreateCheckersDrawingBrush()));
            //MatGR.Children.Add(new SpecularMaterial(Brushes.SkyBlue, 100));
            //BAckMatGR.Children.Add(new DiffuseMaterial(Brushes.Green));
            //BAckMatGR.Children.Add(new DiffuseMaterial(TubeGeometry.CreateCheckersDrawingBrush()));
            
            mGeometry = new GeometryModel3D(mesh, MatGR);
            //mGeometry.BackMaterial = BAckMatGR;
            mGeometry.BackMaterial = MatGR;
            //mGeometry = new GeometryModel3D(mesh, new DiffuseMaterial(Brushes.YellowGreen));
            Transform3DGroup MyTransform3DGroup = new Transform3DGroup();
   

            RotateTransform3D myRotateTransform3D = new RotateTransform3D();
            //rotationY = new AxisAngleRotation3D();
            rotationY.Axis = TubeGeometry.TheAxis;
            //rotationY.Axis = TubeGeometry.GroundAxis;
            rotationY.Angle = 0;
            myRotateTransform3D.Rotation = rotationY;

            MyTransform3DGroup.Children.Add(myRotateTransform3D);
             myMatrixTransform3D = new MatrixTransform3D();
            MyTransform3DGroup.Children.Add(myMatrixTransform3D);

            mGeometry.Transform = MyTransform3DGroup;
            group.Children.Add(mGeometry);


        }

        private void RotateTube(object sender, RoutedEventArgs e)
        {
            Transform3DGroup Tr3D = new Transform3DGroup();

            Transform3DGroup current = (Transform3DGroup)mGeometry.Transform;
            //mGeometry.Transform
            ScaleTransform3D myScaleTransform3D = new ScaleTransform3D();
            myScaleTransform3D.ScaleX = 1;
            myScaleTransform3D.ScaleY = 1;
            myScaleTransform3D.ScaleZ = 1;
            Tr3D.Children.Add(myScaleTransform3D);
           
            //return;
            //for (double  counter = 0; counter < 100; counter = counter + 3)
            //{
                
                
            //    //Tr3D.Children.Add(myScaleTransform3D);


            //    //RotateTransform3D myRotateTransform3D = new RotateTransform3D();
            //    //AxisAngleRotation3D myAxisAngleRotation3d = new AxisAngleRotation3D();
            //    //myAxisAngleRotation3d.Axis = TubeGeometry.TheAxis;
            //    //myAxisAngleRotation3d.Angle = counter;
            //    //myRotateTransform3D.Rotation = myAxisAngleRotation3d;

            //    //Tr3D.Children.Add(myRotateTransform3D);
            //    //ScaleTransform3D myScaleTransform3D = new ScaleTransform3D();
            //    ((ScaleTransform3D)Tr3D.Children[0]).ScaleX = 1 - 1 / counter;
            //    ((ScaleTransform3D)Tr3D.Children[0]).ScaleY = 1 - 1 / counter;
            //    ((ScaleTransform3D)Tr3D.Children[0]).ScaleZ = 1 - 1 / counter;

            //    //Tr3D.Children.Add(myScaleTransform3D);

            //    //mGeometry.Transform = Tr3D;
            //    //mGeometry.Transform.TryTransform();
                
            //}
            RotateTransform3D myRotateTransform3D = new RotateTransform3D();
            AxisAngleRotation3D myAxisAngleRotation3d = new AxisAngleRotation3D();
            myAxisAngleRotation3d.Axis = TubeGeometry.GroundAxis;
            myAxisAngleRotation3d.Angle = 10;
            myRotateTransform3D.Rotation = myAxisAngleRotation3d;

            Tr3D.Children.Add(myRotateTransform3D);
            mGeometry.Transform = Tr3D;
            Vector3D NewPerPendicularAxis = TubeGeometry.TheAxis;
            
        }

        private void Window_MouseEnter_1(object sender, MouseEventArgs e)
        {
            int z = 0;
            UpdateLayout();
        }

    }
    
}
