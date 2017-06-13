using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlowRegime
{
    /// <summary>
    /// Стационарная двухфазная гидродинамическая модель. Рассчитывает объемную долю и градиент давления для заданного сечения (и скорости флюидов).
    /// </summary>

    interface IHydrodynamicModel
    {
        void SetLinerWall(double innerDiameter, double relativeRoughness);
        void SetInclination(double angle);
        void SetHeavyFluid(double density, double viscosity);
        void SetLightFluid(double density, double viscosity);
        void SetSuperficialVelocityies(double light, double heavy);
        void SetInterphaseTension(double tension);

        void Calculate();

        double GetHoldupLight();
        double GetHoldupHeavy();
        double GetPressureGradient();
        double GetVelocityLight();
        double GetVelocityHeavy();
    }
}
