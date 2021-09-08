// self.addEventListener('fetch', event=> {
//     let response=fetch(event.request);
//     event.respondWith(response);
// });
// Plot
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


function calcGeometricFactor(Rl,Rt)
{
    let e = Math.sqrt(1-Math.pow((Rt/Rl),2))
    let L = (1-Math.pow(e,2))/Math.pow(e,2)*((1/(2*e))*Math.log((1+e)/(1-e))-1)
    return L
}
function calcAbsMieGans(em,e1,e2,xi1,xi2,Rl,Rt)
{
    let abs=0;
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
            abs += numerator/denominator;
        }   
    }
    
    return abs;
}

function generatePlot(em,xi1,xi2,Rl,Rt)
{
    myChart.destroy();
    // Initialize the variables
    let lambda = [];
    let absData=[];
    //let lambdaPlotRange = [];
    // let em = 1.76;
    // let Rl = 15;
    // let Rt = 5;
    for(let i=0;i<lda0.length;i=i+1)
    {
        lambda.push(lda0[i]);
    }
    //lambdaPlotRange=[];
    //abs = [];
    for(let j=0;j<lambda.length;j++)
    {
        if(lambda[j]>=400 && lambda[j]<=800)
        {
            //lambdaPlotRange.push(lambda[j]);
            let e1 =eRealAu[j]; //Math.pow(nAu[j],2)-Math.pow(kAu[j],2);
            let e2 =eImagAu[j];// 2*nAu[j]*kAu[j];
            // let xi1=1;
            // let xi2=1;
            absData.push(calcAbsMieGans(em,e1,e2,xi1,xi2,Rl,Rt));
        }
    }
    myChart.config.data.datasets[0].data=absData;
    // myChart.update();
    myChart = new Chart(
        document.getElementById('absorption'),
        config
    );
}

generatePlot(1.76,1,1,15,5);


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
        em = document.getElementById("em").value;
        emScaled= em*100;
        document.getElementById('em_slider').value=emScaled;
    }
    console.log(rTransverse);
    console.log(rLongitude);
    generatePlot(em,1,1,rLongitude,rTransverse); 
}
