let absPlotData=[];
let lambdaPlotRange=[];
const data = {
    labels: lda0,
    datasets: [{
        label: 'Absorption',
        backgroundColor: 'rgb(255, 99, 132)',
        borderColor: 'rgb(255, 99, 132)',
        pointRadius: 1,
        data: absPlotData
    }],
};
const config = {
    type: 'line',
    data,
    options: {}
};

let myChart = new Chart(
    document.getElementById('absorption'),
    config
);

function calcChiralityReal(lambda)
{
    let lambdaC0=document.getElementById('lambdac0').value;
    let gammac=document.getElementById('gammac').value;
    let tauc=document.getElementById('tauc').value;
    let betac=document.getElementById('betac').value;
    let epsilonc0=document.getElementById('epsilonc0').value;

    const CCONST = 299792458; // speed of light
    const HBAR = 1.0545718e-34; // Reduced planck's constant  
    let omega0 = 2*Math.PI*CCONST/lambdaC0;
    let omega = CCONST*2*Math.PI*1e9/lambda;

    let xireal = betac*(2*omega*HBAR* (tauc**2 + (omega-omega0)*(omega+omega0)*HBAR**2))/((tauc**2+((omega-omega)**2)* HBAR**2)*(tauc**2 +((omega+omega0)**2)*HBAR**2))
    return xireal
}

function calcChiralityImag(lambda)
{
    let lambdaC0=document.getElementById('lambdac0').value;
    let gammac=document.getElementById('gammac').value;
    let tauc=document.getElementById('tauc').value;
    let betac=document.getElementById('betac').value;
    let epsilonc0=document.getElementById('epsilonc0').value;

    const CCONST = 299792458; // speed of light
    const HBAR = 1.0545718e-34; // Reduced planck's constant  
    let omega0 = 2*Math.PI*CCONST/lambdaC0;
    let omega = CCONST*2*Math.PI*1e9/lambda;

    let xiimag = betac*tauc * (-( 1/(tauc**2 + ((omega + omega0)**2)*HBAR**2))-1/(tauc**2 + (omega*HBAR - omega0*HBAR)**2));

    return xiimag
}
function calcGeometricFactor(Rl,Rt)
{
    let e = Math.sqrt(1-Math.pow((Rt/Rl),2))
    let L = (1-Math.pow(e,2))/Math.pow(e,2)*((1/(2*e))*Math.log((1+e)/(1-e))-1)
    return L
}

function calcAbsMieGans(em,e1,e2,xi1,xi2,Rl,Rt)
{
    let absL=0;
    let Lx = calcGeometricFactor(Rl,Rt);
    let Ly=(1-Lx)/2;
    let Lz=Ly;
    let L =[Lx,Ly,Lz];
    for(let i=0;i<3;i++)
    {
        let ps=L[i];
        let fs=(em+(e1-em)*ps)*xi2+ps*e2*xi1;
        let cs = [1,1]
        let ct = [1,-1]
        for(let j=0;j<2;j++)
        {
            let numerator=Math.pow(cs[j],2)*Math.sqrt(em)*e2+2*cs[j]*ct[j]*fs;
            let denominator=Math.pow(((ps-1)*em-ps*e1),2)+Math.pow((ps*e2),2);
            absL += numerator/denominator;
        }   
    }
    return absL;
}

function generatePlot(em,Rl,Rt,metal)
{
    myChart.destroy();
    let e1=0;
    let e2=0;
    let lambda = [];
    let absData=[];

    

    for(let i=0;i<lda0.length;i=i+1)
    {
        lambda.push(lda0[i]);
    }
    for(let j=0;j<lambda.length;j++)
    {
        if(lambda[j]>=400 && lambda[j]<=1000)
        {
            if(metal =='Au')
            {
                e1 =eRealAu[j]; 
                e2 =eImagAu[j];
            }
            else if (metal =='Ag')
            {
                e1 =eRealAg[j]; 
                e2 =eImagAg[j];
            }
            let xi1=calcChiralityReal(lambda[j]);
            let xi2=calcChiralityImag(lambda[j]);

            //effective medium parameters should be pushed.
            absData.push(calcAbsMieGans(em,e1,e2,xi1,xi2,Rl,Rt));
        }
    }
    myChart.config.data.datasets[0].data=absData;
    myChart = new Chart(
        document.getElementById('absorption'),
        config
    );
}

generatePlot(1.76,15,5,'Au');


function updatePlot(slider)
{
    let emScaled = 0;
    let em = 0;
    let rTransverse = document.getElementById('Rt').value;
    let rLongitude = document.getElementById('Rl').value;
    if (slider)
    {
        emScaled = document.getElementById('em_slider').value;
        em= emScaled/100;
        document.getElementById("em").value=em;
    }
    else
    {
        em = Number(document.getElementById("em").value);
        emScaled= em*100;
        document.getElementById('em_slider').value=emScaled;
    }


    let metal=document.querySelector('input[name="metal"]:checked').value;
    console.log(metal);
    generatePlot(em,rLongitude,rTransverse,metal); 
}