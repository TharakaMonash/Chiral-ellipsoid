let absLeftPlotData = [];
let absRightPlotData = [];
let cdPlotData = [];
let lambdaPlotRange = [];
const dataAbs = {
    labels: lda0,
    datasets: [{
        label: 'Absorption Left',
        backgroundColor: 'rgb(255, 99, 132)',
        borderColor: 'rgb(255, 99, 132)',
        pointRadius: 1,
        data: absLeftPlotData
    },
    {
        label: 'Absorption Right',
        backgroundColor: 'rgb(0, 0, 0)',
        borderColor: 'rgb(0, 0, 0)',
        pointRadius: 1,
        data: absRightPlotData
    }
    ],
};
const dataCd = {
    labels: lda0,
    datasets: [
        {
            label: 'CD',
            backgroundColor: 'rgb(0, 0, 0)',
            borderColor: 'rgb(0, 0, 0)',
            pointRadius: 1,
            data: cdPlotData
        }
    ],
};
const configAbs = {
    type: 'line',
    data: dataAbs,
    options: {
        scales: {
            x: {
                title: {
                    display: true,
                    text: 'Wavelength'
                }
            },
            y: {
                title: {
                    display: true,
                    text: 'Absorption'
                }
            }
        }
    }
};

const configCd = {
    type: 'line',
    data: dataCd,
    options: {
        scales: {
            x: {
                title: {
                    display: true,
                    text: 'Wavelength'
                }
            },
            y: {
                title: {
                    display: true,
                    text: 'CD'
                }
            }
        }
    }
};

let myChart = new Chart(
    document.getElementById('absorption'),
    configAbs
);

let cdChart = new Chart(
    document.getElementById('cd'),
    configCd
);

function calcChirality(lambda) {
    let lambdaC0 = Number(document.getElementById('lambdac0').value);
    let gammac = Number(document.getElementById('gammac').value);
    let tauc = Number(document.getElementById('tauc').value);
    let betac = Number(document.getElementById('betac').value);
    let ecm = Number(document.getElementById('epsilonc0').value);

    const CCONST = 299792458; // speed of light
    const HBAR = 1.0545718e-34; // Reduced planck's constant  
    let omega0 = 2 * Math.PI * CCONST / lambdaC0;
    let omega = CCONST * 2 * Math.PI * 1e9 / lambda;

    let xireal = betac * (2 * omega * HBAR * (tauc ** 2 + (omega - omega0) * (omega + omega0) * (HBAR ** 2))) / ((tauc ** 2 + ((omega - omega0) ** 2) * HBAR ** 2) * (tauc ** 2 + ((omega + omega0) ** 2) * HBAR ** 2));

    let xiimag = betac * tauc * (-(1 / (tauc ** 2 + ((omega + omega0) ** 2) * HBAR ** 2)) - 1 / (tauc ** 2 + (omega * HBAR - omega0 * HBAR) ** 2));

    let eCreal =ecm+(2 *gammac *HBAR *omega0* (HBAR**2 *(omega0**2-omega**2)+tauc**2))/((HBAR**2 *(omega-omega0)**2+tauc**2)* (HBAR**2 *(omega+omega0)**2+tauc**2));

    let eCimag = gammac *tauc* (-1*(1/(HBAR**2* (omega+omega0)**2+tauc**2))+1/((HBAR *omega-HBAR *omega0)**2+tauc**2))
    let chiral={
        eC:{
            real:eCreal,
            imag:eCimag
        },
        xi : {
            real: xireal,
            imag: xiimag
        }
    }
    return chiral;
}

function calcGeometricFactor(Rl, Rt) {
    let e = Math.sqrt(1 - Math.pow((Rt / Rl), 2))
    let L = (1 - Math.pow(e, 2)) / Math.pow(e, 2) * ((1 / (2 * e)) * Math.log((1 + e) / (1 - e)) - 1)
    return L
}


function calcAbsN(eMedium, lambda, e1, e2, longRadius, transRadius, xi1, xi2, N) {
    // Absorption of a longitudinal ellipsoid with left polarized light   
    let AbsLlong = calcAbsOneEllipsoid(eMedium, lambda, e1, e2, longRadius, transRadius, xi1, xi2, -1, 1 / Math.sqrt(2), 1);
    // Absorption of a longitudinal ellipsoid with right polarized light  
    let AbsRlong = calcAbsOneEllipsoid(eMedium, lambda, e1, e2, longRadius, transRadius, xi1, xi2, 1, 1 / Math.sqrt(2), 1);
    // Absorption of a transverse ellipsoid with left polarized light
    let AbsLtrans = calcAbsOneEllipsoid(eMedium, lambda, e1, e2, longRadius, transRadius, xi1, xi2, -1, 1 / Math.sqrt(2), 0);
    // Absorption of a transverse ellipsoid with right polarized light
    let AbsRtrans = calcAbsOneEllipsoid(eMedium, lambda, e1, e2, longRadius, transRadius, xi1, xi2, 1, 1 / Math.sqrt(2), 0);
    // Absorption of N ellipsoids with left polarized light   
    let AbsL = N * (AbsLlong + AbsLlong + AbsLtrans) / 3;
    // Absorption of N ellipsoids with right polarized light   
    let AbsR = N * (AbsRlong + AbsRlong + AbsRtrans) / 3;
    let cd = AbsR - AbsL;
    let AbsObject =
    {
        left: AbsL,
        right: AbsR,
        cd: cd
    };
    return AbsObject;
}


function calcAbsOneEllipsoid(eMedium, lambda, e1, e2, longRadius, transRadius, xi1, xi2, C, K, orientation) {
    // Geometric factor calculation
    let R = transRadius / longRadius;
    let e = Math.sqrt(1 - R ** 2);
    let L = ((1 - e ** 2) / e ** 2) * ((1 / (2 * e)) * Math.log((1 + e) / (1 - e)) - 1);
    if (orientation == 1) // "long"
    {
        P = [L, (1 - L) / 2, (1 - L) / 2];
    }
    else if (orientation == 0) {
        P = [(1 - L) / 2, (1 - L) / 2, L];
    }

    let N = 1;
    let SL = 299792458; // Speed of light
    let V = 4 / 3 * Math.PI * longRadius * transRadius * transRadius * (1e-9) ** 3;

    let temp = 0;
    for (let i = 0; i < 2; i++) {
        temp = temp + (Math.sqrt(eMedium) * e2 - 2 * C * ((eMedium + P[i] * (e1 - eMedium)) * xi2 + P[i] * e2 * xi1)) / (((P[i] - 1) * eMedium - P[i] * e1) ** 2 + (P[i] * e2) ** 2);
    }
    abs = N * V * (K ** 2) * SL * (eMedium ** (3 / 2)) * temp / lambda;
    return abs;
}

function effectiveMedium(e1, e2, eC1,eC2, frac) {
    
    let effReal = (9 *e2**2 *eC1* frac-e1 *(-1+frac) *(eC1**2+(2 *e2+eC2)**2+2* (eC1**2+(e2-eC2)**2) *frac)-2 *e1**3 *(-2+frac+frac**2)+e1**2 *eC1 *(4+frac+4* frac**2))/((eC1**2+eC2**2) *(-1+frac)**2+e1**2* (2+frac)**2+e2**2 *(2+frac)**2-2* e1 *eC1* (-2+frac+frac**2)-2 *e2 *eC2* (-2+frac+frac**2));
    let effImag = (4* e1 *e2* eC1 *(-1+frac)**2+e1**2 *(9 *eC2 *frac-2* e2 *(-2+frac+frac**2))+e2* (-((eC1**2+eC2**2) *(-1+frac)* (1+2* frac))-2 *e2**2 *(-2+frac+frac**2)+e2 *eC2 *(4+frac+4 *frac**2)))/((eC1**2+eC2**2) *(-1+frac)**2+e1**2*(2+frac)**2+e2**2* (2+frac)**2-2 *e1* eC1* (-2+frac+frac**2)-2* e2 *eC2* (-2+frac+frac**2));
    let eeff = {
        real: effReal,
        imag: effImag
    };
    return eeff;
}

function effectiveChirality(xi1, xi2, e1,e2, eC1,eC2, frac) {
    
    let xiEffReal = (3 *frac* (e1**2 *(2+frac) *xi1+e2* (eC2-eC2* frac+e2 *(2+frac)) *xi1+e2 *eC1 *(-1+frac) *xi2-e1* (-1+frac) *(eC1* xi1+eC2* xi2)))/((eC1**2+eC2**2)* (-1+frac)**2+e1**2 *(2+frac)**2+e2**2 *(2+frac)**2-2* e1* eC1 *(-2+frac+frac**2)-2* e2* eC2 *(-2+frac+frac**2)) ;
    let xiEffImag = (3* frac *(e2**2 *(2+frac)* xi2-e2* (-1+frac) *(eC1 *xi1+eC2 *xi2)+e1 *(e1* (2+frac) *xi2+(-1+frac)* (eC2* xi1-eC1* xi2))))/((eC1**2+eC2**2)* (-1+frac)**2+e1**2 *(2+frac)**2+e2**2 *(2+frac)**2-2 *e1 *eC1 *(-2+frac+frac**2)-2* e2* eC2 *(-2+frac+frac**2));
    let xiEff = {
        real: xiEffReal,
        imag: xiEffImag
    };
    return xiEff;
}

function generatePlot(em, Rl, Rt, metal) {
    myChart.destroy();
    cdChart.destroy();
    let e1 = 0;
    let e2 = 0;
    let lambda = [];
    let absLeftData = [];
    let absRightData = [];
    let cdData = [];


    for (let i = 0; i < lda0.length; i = i + 1) {
        lambda.push(lda0[i]);
    }
    for (let j = 0; j < lambda.length; j++) {
        if (lambda[j] >= 400 && lambda[j] <= 1000) {
            if (metal == 'Au') {
                e1 = eRealAu[j];
                e2 = eImagAu[j];
            }
            else if (metal == 'Ag') {
                e1 = eRealAg[j];
                e2 = eImagAg[j];
            }
            let chiral = calcChirality(lambda[j]);
            let xi1 = chiral.xi.real;
            let xi2 = chiral.xi.imag;
            let eC1 = chiral.eC.real;
            let eC2 = chiral.eC.imag;
            let N = 1e12;
            let frac =Number(document.getElementById('frac').value);
            //effective medium parameters should be pushed.
            let eEff = effectiveMedium(e1,e2,eC1,eC2,frac);
            let xiEff = effectiveChirality(xi1,xi2,e1,e2,eC1,eC2,frac);
            let eEffR = eEff.real;
            let eEffI = eEff.imag;
            let xiEff1 = xiEff.real;
            let xiEff2 = xiEff.imag;
            //absData.push(calcAbsMieGans(em,e1,e2,xi1,xi2,Rl,Rt));
            let object = calcAbsN(em, lambda[j], eEffR, eEffI, Rl, Rt, xiEff1, xiEff2, N);
            absLeftData.push(object.left);
            absRightData.push(object.right);
            cdData.push(object.cd);
        }
    }
    myChart.config.data.datasets[0].data = absLeftData;
    myChart.config.data.datasets[1].data = absRightData;
    cdChart.config.data.datasets[0].data = cdData;
    myChart = new Chart(
        document.getElementById('absorption'),
        configAbs
    );
    cdChart = new Chart(
        document.getElementById('cd'),
        configCd
    );
}

function updatePlot(slider) {
    let emScaled = 0;
    let em = 0;
    let rTransverse = document.getElementById('Rt').value;
    let rLongitude = document.getElementById('Rl').value;
    if (slider) {
        emScaled = document.getElementById('em_slider').value;
        em = emScaled / 100;
        document.getElementById("em").value = em;
    }
    else {
        em = Number(document.getElementById("em").value);
        emScaled = em * 100;
        document.getElementById('em_slider').value = emScaled;
    }

    let metal = document.querySelector('input[name="metal"]:checked').value;
    console.log(metal);
    generatePlot(em, rLongitude, rTransverse, metal);
}

generatePlot(1.76, 15, 5, 'Au');