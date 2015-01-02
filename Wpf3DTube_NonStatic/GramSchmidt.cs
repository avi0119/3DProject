using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Media.Media3D;
using System.Collections.Generic;

namespace Wpf3DTube_NonStatic
{
    public class GramSchmidt
    {
        private List<Vector3D> OriginalBasis = new List<Vector3D>();
        private List<Vector3D> OrthogonalBasis = new List<Vector3D>();
        public List<Vector3D> OrthonormalBasis 
        { 
            get
            {
                return OrthogonalBasis;
            }

        }
        public GramSchmidt()
        { 
        }
        public GramSchmidt(bool normalizeBasis, params Vector3D[] list)
        {
            for (int i = 0; i < list.Length; i++)
            {

                OriginalBasis.Add(list[i]);
            }
            createOrthogonalBasis( normalizeBasis);

        }
        private void createOrthogonalBasis(bool normalizeBasis)
        {
            Vector3D currVector;
            Vector3D currVectorMinusProjection;
            Vector3D currentProjection;
            for (int i = 0; i < OriginalBasis.Count; i++)
            {

                currVector=OriginalBasis[i];
                if (i == 0)
                {
                    OrthogonalBasis.Add(currVector);
                }
                else
                {
                    currentProjection = returnProjetionOfVectorPerGivenBasis(currVector, OrthogonalBasis);
                    currVectorMinusProjection = currVector - currentProjection;
                    OrthogonalBasis.Add(currVectorMinusProjection);
                }
            }
            if (normalizeBasis == true)
            {
                for (int i = 0; i < OrthogonalBasis.Count; i++)
                {
                    Vector3D v = OrthogonalBasis[i];
                    v.Normalize();
                    OrthogonalBasis[i] = v;

                }
            }
        }
        private Vector3D returnProjetionOfVectorPerGivenBasis(Vector3D v, List<Vector3D> anOrthogonalBasis)
        {
            Vector3D cumulative=new Vector3D();
            Vector3D currentBasisVector;
            Vector3D currentProjection;
            for (int i = 0; i < anOrthogonalBasis.Count; i++)
            {
                currentBasisVector = OrthogonalBasis[i];
                currentProjection = (Vector3D.DotProduct(v, currentBasisVector)/(currentBasisVector.Length*currentBasisVector.Length))*currentBasisVector;
                cumulative=cumulative+currentProjection;
            }
            return cumulative;
        }
    }
}
