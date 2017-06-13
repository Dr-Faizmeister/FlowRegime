using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace FlowRegime
{
    public partial class Form1 : Form
    {
        public byte FlowRegime;

        public float D, theta, rho_l, rho_g, visc_l, visc_g, Q_l, Q_g, sigma, h_pr;
        public double V_sg, V_sl;

        public Form1()
        {
            InitializeComponent();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            D = float.Parse(textBox1.Text);
            rho_l = float.Parse(textBox2.Text);
            rho_g = float.Parse(textBox3.Text);
            visc_l = float.Parse(textBox4.Text);
            visc_g = float.Parse(textBox5.Text);
            theta = float.Parse(textBox6.Text);
            sigma = float.Parse(textBox7.Text);
            h_pr = float.Parse(textBox8.Text);
            Q_l = float.Parse(textBox9.Text);
            Q_g = float.Parse(textBox10.Text);

            V_sl = 4 * Q_l / (Math.PI * D * D * 86400);
            V_sg = 4 * Q_g / (Math.PI * D * D * 86400);

            FlowPattern flowPattern = new FlowPattern(D, rho_l, rho_g, visc_l, visc_g, sigma, h_pr);

            FlowRegime = flowPattern.GetFlowRegime(V_sl, V_sg, theta);

            string[] flowRegimeName =
        {
            "N/A",
            "DispersedBubble",
            "StratifiedSmooth",
            "StratifiedWavy",
            "AnnularMist",
            "Bubble",
            "Slug",
            "ElongatedBubble",
            "Froth"
        };
            textBox13.Text = flowRegimeName[FlowRegime];
            textBox11.Text = (flowPattern.GetPressureGradient(V_sl, V_sg, theta) / 101325).ToString();
            textBox12.Text = (flowPattern.GetLiquidHoldup(V_sl, V_sg, theta)).ToString();
        }
    }
}
