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
    public class BASIS
    {
        public BASIS()
        {

        }
        static public  Vector3D[] findOrthonormalBasisToVector (Vector3D n)
        {
            Point3D p1 = new Point3D();
            #region findingPointOnPlane
            double res;
            Point3D somePointInPlane=new Point3D();
            Plane nPlane=new Plane(n, p1);
            PlaneEquation equation = nPlane.Equation;

            EquationParameter EP = new EquationParameter();

            EP.X = 1.0;
            EP.Z = 0.0;
           
            if (nPlane.Equation.solveForRemainingVariable(EP, out res) == true)
            {
                // solving for y
                somePointInPlane = new Point3D((double)EP.X, res, (double)EP.Z);
            }
            else
            {
                EP.RotatePrameters();
                if (nPlane.Equation.solveForRemainingVariable(EP, out res) == true)
                {
                    // solving fo z
                    somePointInPlane = new Point3D((double)EP.X, (double)EP.Y, res);
                }
                else
                {
                    EP.RotatePrameters();
                    if (nPlane.Equation.solveForRemainingVariable(EP, out res) == true)
                    {
                        // solving fo x
                        somePointInPlane = new Point3D(res, (double)EP.Y, (double)EP.Z);
                    }
                }
            }
            #endregion
            Vector3D v1=p1-somePointInPlane;
            Vector3D v2=Vector3D.CrossProduct(v1,n);
            v1.Normalize();
            v2.Normalize();
            Vector3D[] basis=new Vector3D[]{v1,v2};
            //GramSchmidt GM=new GramSchmidt (true,v);
            //Vector3D[] basis= GM.OrthonormalBasis.ToArray();
            return basis;
        }
    }
}
