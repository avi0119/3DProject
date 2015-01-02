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
using System.ComponentModel;


namespace Wpf3DTube_NonStatic
{
    class Point3DInfo
    {
        Point3D prPoint;
        double prcolorIntensity;
        DateTime prTimeStamp;
        SolidColorBrush prslb;
        MyColor prMyColor;
        List<GeometryModel3D> prlGeoMo3D;

        public List<GeometryModel3D> lGeoMo3D
        {
            get
            {
                return prlGeoMo3D;
            }
            set
            {
                prlGeoMo3D = value;
            }
        }

        Material prMaterial;
        public Material ShapesMaterial
        {
            get
            {
                return prMaterial;
            }
            set
            {
                prMaterial = value;
            }
        }



        MyMaterial prMyMaterial;
        public MyMaterial mymaterial
        {
            get
            {
                return prMyMaterial;
            }
            set
            {
                prMyMaterial = value;
            }
        }

        public MyColor myColor
        {
            get
            {
                return prMyColor;
            }
            set
            {
                prMyColor = value;
            }
        }

        public SolidColorBrush SLB
        {
            get
            {
                return prslb;
            }
            set
            {
                prslb = value;
            }
        }
        public DateTime TimeStamp
        {
            get
            {
                return prTimeStamp;
            }
            set
            {
                prTimeStamp = value;
            }
        }
        public double ColorIntensity
        {
            get
            {
                return prcolorIntensity;
            }
            set
            {
                prcolorIntensity = value;
            }
        }
        public Point3D Point
        {
            get
            {
                return prPoint;
            }
            set
            {
                prPoint = value;
            }
        }
        public Point3DInfo()
        {

        }
        public Point3DInfo(Point3D aPoint)
        {
            Point = aPoint;
            TimeStamp = DateTime.Now;
            //ColorIntensity = aIntensity;
        }

    }
    class Trail3DPoints
    {
        const double powerbase = 0.93;
        Model3DGroup prmyModel3DGroup;
        int prPointMeter=0;
        List<Point3DInfo> prTrail = new List<Point3DInfo>();
        double prThreshold;
        int prFrequencyOfRecording;
        public double Threshold
        {
            get
            {
                return prThreshold;
            }
            set
            {
                prThreshold = value;
            }
        }
        public Model3DGroup myModel3DGroup
        {
            get
            {
                return prmyModel3DGroup;
            }
            set
            {
                prmyModel3DGroup = value;
            }
        }
        public int FrequencyOfRecording
        {
            get
            {
                return prFrequencyOfRecording;
            }
            set
            {
                prFrequencyOfRecording = value;
            }
        }
        public int PointMeter
        {
            get
            {
                return prPointMeter;
            }
            set
            {
                prPointMeter = value;
            }
        }
        public Trail3DPoints()
        {

        }
        public Trail3DPoints(int PointsToRecordFrequency,double thresholdBelowWhichToRemove,Model3DGroup model3dgroup)
        {
            Threshold = thresholdBelowWhichToRemove;
            FrequencyOfRecording = PointsToRecordFrequency;
            myModel3DGroup = model3dgroup;
        }
        private bool areTwoPointEqual(Point3D A,Point3D B)
        {
            return ( (A.X==B.X ) && ( A.Y==B.Y ) && ( A.Z==B.Z ));
        }
        public void AddPointIfApplicable(Point3D aActualPoint)
        {
            //bool isdiff=false;
            Point3DInfo aPointInfo = new Point3DInfo(aActualPoint);
            aPointInfo.ColorIntensity = 1;
            PointMeter = PointMeter % FrequencyOfRecording;
            PointMeter = PointMeter + 1;
            //Point3D lastpoint;
            
            if ((PointMeter == 1) )
            {
                //if (prTrail.Count > 0)
                //{
                //    lastpoint = prTrail[prTrail.Count-1].Point;
                //    isdiff = areTwoPointEqual(aActualPoint, lastpoint);
                //}
                //else
                //{
                //   isdiff= false;
                //}
                //if (isdiff == false)
                //{
                      List<GeometryModel3D> lGeoMo3D = new List<GeometryModel3D>();
                      aPointInfo.lGeoMo3D = lGeoMo3D;
                     SolidColorBrush slb = new SolidColorBrush(Colors.Blue);
                     aPointInfo.SLB = slb;
                     //MyColor a = bindBrushColorPropertyToColor(slb);
                     //aPointInfo.myColor = a;
                     Material mtOfAllSidesOFcube= new DiffuseMaterial((new SolidColorBrush(Colors.Blue)));
                     aPointInfo.ShapesMaterial = mtOfAllSidesOFcube;
                     MyMaterial myNewmaterial = new MyMaterial(new DiffuseMaterial(new SolidColorBrush(Colors.Blue)));
                     aPointInfo.mymaterial = myNewmaterial;
                    prTrail.Add(aPointInfo);
                    CreateCubeWithinGivenModel3DGroup(lGeoMo3D,mtOfAllSidesOFcube, myNewmaterial, myModel3DGroup, aActualPoint.X, aActualPoint.Y, aActualPoint.Z, slb, 0.75);
                //}
                //else
                //{
                //    isdiff = false;
                //    bool h = isdiff;
                //}
               

            }
            RemoveColorlessPoints();
            AdjustTrail();

        }
        public void AdjustTrail()
        {
            int count = prTrail.Count;
            int pointcounter = count;
            int internalcounter = 0;
            
            double currentIntensity = 0;
            foreach (Point3DInfo a in prTrail)
            {
                internalcounter = internalcounter + 1;
                DateTime now=DateTime.Now;
                TimeSpan timeElapsed = now - a.TimeStamp;
                int NumberOfMillisecomdsElapsed = timeElapsed.Seconds * 1000 + timeElapsed.Milliseconds * 1;
                double raisetopower = ((double)NumberOfMillisecomdsElapsed) / 100;
                //currentIntensity = Math.Pow(powerbase, raisetopower);
                currentIntensity = Math.Pow(powerbase, (count-internalcounter)*.35);
                a.ColorIntensity = currentIntensity;
                //a.myColor.color = createBlueColorWihhGivenAlphaChannel(currentIntensity);
                //a.mymaterial.material = createMateruakWithBlueColorWihhGivenAlphaChannel(currentIntensity);
                //a.ShapesMaterial.
                Material nm = createMateruakWithBlueColorWihhGivenAlphaChannel(currentIntensity);
                foreach (GeometryModel3D item in a.lGeoMo3D)
                {
                    item.Material =nm;
                    item.BackMaterial = nm;
                }
                pointcounter = pointcounter - 1;
                //a.SLB.Color.B = 255;
            }
                
        }
        public void RemoveColorlessPoints_old2()
        {
            int count;
            int counter;
            double currentIntensity = 0;
            //Point3DInfo a;
            List<int> tobremoved = new List<int>();
            foreach (Point3DInfo a in prTrail)
            {
                currentIntensity = a.ColorIntensity;
                if ((currentIntensity < Threshold))
                {
                    prTrail.Remove(a);
                    RemoveOldestPointCreated();


                }
            }



            //count = prTrail.Count;
            //for (counter = 0; counter < count; counter = counter + 1)
            //{
            //    a = prTrail[counter];
            //    currentIntensity = a.ColorIntensity;
            //    if ((currentIntensity < Threshold))
            //    {
            //        prTrail.RemoveAt(counter);
            //        RemoveOldestPointCreated();


            //    }
            //}
            //if (prTrail.Count > 100)
            //{
            //    RemoveLastPoint(true);
            //    prTrail.Clear();
            //}
        }
        public void RemoveColorlessPoints ()
        {
            int count ;
            int counter;
            double currentIntensity = 0;
            Point3DInfo a ;
            count = prTrail.Count;
            List<int> tobremoved = new List<int>();
            for (counter = 0; counter < count; counter = counter + 1)
            {
                a=prTrail[counter];
                currentIntensity = a.ColorIntensity;
                if ((currentIntensity < Threshold) )
                {
                    tobremoved.Add(counter);
                    //prTrail.RemoveAt(counter);
                    //RemoveOldestPointCreated();


                }
            }
            for (counter = tobremoved.Count; counter >0; counter = counter - 1)
            {

                int indetoremove = tobremoved[counter-1];
                //a = prTrail[indetoremove];


                    prTrail.RemoveAt(indetoremove);
                    RemoveOldestPointCreated();


                
            }
            //if (prTrail.Count > 100)
            //{
            //    RemoveLastPoint(true);
            //    prTrail.Clear();
            //}
        }
        protected static Vector3D CalculateNormal(Point3D p0, Point3D p1, Point3D p2)
        {
            Vector3D v0 = new Vector3D(p1.X - p0.X, p1.Y - p0.Y, p1.Z - p0.Z);
            Vector3D v1 = new Vector3D(p1.X - p2.X, p1.Y - p2.Y, p2.Z - p1.Z);

            return Vector3D.CrossProduct(v0, v1);
        }
        protected Model3DGroup CreateTriangle(Point3D p0, Point3D p1, Point3D p2, SolidColorBrush slb, MyMaterial mmat, Material mtOfAllSidesOFcube,List<GeometryModel3D> lGeoMo3D)
        {
            MeshGeometry3D mesh = new MeshGeometry3D();
           
            mesh.Positions.Add(p0);
            mesh.Positions.Add(p1);
            mesh.Positions.Add(p2);
            mesh.TriangleIndices.Add(0);
            mesh.TriangleIndices.Add(1);
            mesh.TriangleIndices.Add(2);

            Vector3D normal = CalculateNormal(p0, p1, p2);
            mesh.Normals.Add(normal);
            mesh.Normals.Add(normal);
            mesh.Normals.Add(normal);

            //Material material = new DiffuseMaterial(new SolidColorBrush(Colors.Blue));
            Material material = mtOfAllSidesOFcube;// new DiffuseMaterial(new SolidColorBrush(Colors.Blue));
            GeometryModel3D model = new GeometryModel3D(mesh, material);
            lGeoMo3D.Add(model);
            
            Model3DGroup group = new Model3DGroup();
            group.Children.Add(model);
            //bindMaterialPropertyOfModel3dGroupToMaterial(model, mmat);
            return group;
        }
        protected void CreateCubeWithinGivenModel3DGroup(List<GeometryModel3D> lGeoMo3D,Material mtOfAllSidesOFcube, MyMaterial mmat, Model3DGroup existinModel3DGroup, double x, double y, double z,SolidColorBrush slb, double cubesize = 5)
        {
            //Model3DGroup cube = new Model3DGroup();
            //SolidColorBrush slb = new SolidColorBrush(Colors.Blue);
            double A = cubesize;
            double B = cubesize;
            double C = cubesize;

            Point3D p0 = new Point3D(0 + x, 0 + y, 0 + z);
            Point3D p1 = new Point3D(A + x, 0 + y, 0 + z);
            Point3D p2 = new Point3D(A + x, 0 + y, C + z);
            Point3D p3 = new Point3D(0 + x, 0 + y, C + z);
            Point3D p4 = new Point3D(0 + x, B + y, 0 + z);
            Point3D p5 = new Point3D(A + x, B + y, 0 + z);
            Point3D p6 = new Point3D(A + x, B + y, C + z);
            Point3D p7 = new Point3D(0 + x, B + y, C + z);

            //front
            existinModel3DGroup.Children.Add(CreateTriangle(p3, p2, p6, slb, mmat,mtOfAllSidesOFcube,lGeoMo3D));
            existinModel3DGroup.Children.Add(CreateTriangle(p3, p6, p7, slb, mmat,mtOfAllSidesOFcube,lGeoMo3D));

            //right
            existinModel3DGroup.Children.Add(CreateTriangle(p2, p1, p5, slb, mmat,mtOfAllSidesOFcube,lGeoMo3D));
            existinModel3DGroup.Children.Add(CreateTriangle(p2, p5, p6, slb, mmat,mtOfAllSidesOFcube,lGeoMo3D));

            //back
            existinModel3DGroup.Children.Add(CreateTriangle(p1, p0, p4, slb, mmat,mtOfAllSidesOFcube,lGeoMo3D));
            existinModel3DGroup.Children.Add(CreateTriangle(p1, p4, p5, slb, mmat,mtOfAllSidesOFcube,lGeoMo3D));

            //left
            existinModel3DGroup.Children.Add(CreateTriangle(p0, p3, p7, slb, mmat,mtOfAllSidesOFcube,lGeoMo3D));
            existinModel3DGroup.Children.Add(CreateTriangle(p0, p7, p4, slb, mmat,mtOfAllSidesOFcube,lGeoMo3D));

            //top
            existinModel3DGroup.Children.Add(CreateTriangle(p7, p6, p5, slb, mmat,mtOfAllSidesOFcube,lGeoMo3D));
            existinModel3DGroup.Children.Add(CreateTriangle(p7, p5, p4, slb, mmat,mtOfAllSidesOFcube,lGeoMo3D));

            //bottom
            existinModel3DGroup.Children.Add(CreateTriangle(p2, p3, p0, slb, mmat,mtOfAllSidesOFcube,lGeoMo3D));
            existinModel3DGroup.Children.Add(CreateTriangle(p2, p0, p1, slb, mmat,mtOfAllSidesOFcube,lGeoMo3D));

            //ModelVisual3D model = new ModelVisual3D();
            //model.Content = cube;
            //return model;
        }
        protected void RemoveOldestPointCreated (bool removeall=false) 
        {
            bool exit = false;
            int count;
            while (exit == false)
            {
                count = myModel3DGroup.Children.Count;
                if (count > 0)
                {
                    for (int i = 12; i >0; i = i - 1)
                    {
                        //myModel3DGroup.Children[count - 1].
                        myModel3DGroup.Children.RemoveAt(i-1);
                    }
                    //myModel3DGroup.Children.cl
                    if (removeall == false) { exit = true; }
                }
                else
                {
                    exit = true;
                }
            }
        }

        private MyMaterial bindMaterialPropertyOfModel3dGroupToMaterial(GeometryModel3D model, MyMaterial myNewmaterial)
        {
            //Model3DGroup ff = new Model3DGroup();
           
            Binding myBinding = new Binding("TheMaterial");
            myBinding.Source = myNewmaterial;
            BindingOperations.SetBinding(model, GeometryModel3D.MaterialProperty, myBinding);

            return myNewmaterial;
        }

        private MyColor bindBrushColorPropertyToColor(SolidColorBrush slb)
        {
            MyColor myNewColor = new MyColor(new Color());
            Binding myBinding = new Binding("TheColor");
            myBinding.Source = myNewColor;
            BindingOperations.SetBinding(slb, SolidColorBrush.ColorProperty, myBinding);
            return myNewColor;
        }

        private MyColorIntensity bindIntensityOfColorToSoliColorBrush_v2(SolidColorBrush slb)
        {
            MyColorIntensity myDataObject = new MyColorIntensity(1);
            Binding myBinding = new Binding("IntensityAlpha");
            myBinding.Source = myDataObject;
            return myDataObject;
            //BindingOperations.SetBinding(slb, SolidColorBrush., myBinding);
        }
        private System.Windows.Media.Color  createBlueColorWihhGivenAlphaChannel(double Intentisty)
        {
            System.Windows.Media.Color a=new System.Windows.Media.Color();
            a.ScA = (float)Intentisty;
            a.B = 255;
            a.R = 0;
            a.G = 255;
            return a;
        }
        private Material createMateruakWithBlueColorWihhGivenAlphaChannel(double Intentisty)
        {
            Color newcolor = createBlueColorWihhGivenAlphaChannel(Intentisty);
            Material a = new DiffuseMaterial(new SolidColorBrush(newcolor));

            return a;
        }
    }
    public class MyColorIntensity : INotifyPropertyChanged
    {
        private float intensityAlpha;

        public MyColorIntensity() { }

        public MyColorIntensity(float intensityinit)
        {
            intensityAlpha = intensityinit;
        }

        public float IntensityAlpha
        {
            get { return intensityAlpha; }
            set
            {
                intensityAlpha = value;
                OnPropertyChanged("IntensityAlpha");
            }
        }

        public event PropertyChangedEventHandler PropertyChanged;

        private void OnPropertyChanged(string info)
        {
            PropertyChangedEventHandler handler = PropertyChanged;
            if (handler != null)
            {
                handler(this, new PropertyChangedEventArgs(info));
            }
        }
    }

    public class MyColor : INotifyPropertyChanged
    {
        private Color aColor;

        public MyColor() { }

        public MyColor(Color theColor)
        {
            aColor = theColor;
        }

        public Color color
        {
            get { return aColor; }
            set
            {
                aColor = value;
                OnPropertyChanged("TheColor");
            }
        }

        public event PropertyChangedEventHandler PropertyChanged;

        private void OnPropertyChanged(string info)
        {
            PropertyChangedEventHandler handler = PropertyChanged;
            if (handler != null)
            {
                handler(this, new PropertyChangedEventArgs(info));
            }
        }
    }


    public class MyMaterial : INotifyPropertyChanged
    {
        private Material amaterial;

        public MyMaterial() { }

        public MyMaterial(Material thematerial)
        {
            amaterial = thematerial;
        }

        public Material material
        {
            get { return amaterial; }
            set
            {
                amaterial = value;
                OnPropertyChanged("TheMaterial");
            }
        }

        public event PropertyChangedEventHandler PropertyChanged;

        private void OnPropertyChanged(string info)
        {
            PropertyChangedEventHandler handler = PropertyChanged;
            if (handler != null)
            {
                handler(this, new PropertyChangedEventArgs(info));
            }
        }
    }
    
}
