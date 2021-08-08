// self.addEventListener('fetch', event=> {
//     let response=fetch(event.request);
//     event.respondWith(response);
// });

function calcAbs(lda)
{
    return lda**2;
}

// Initialize the variables
let lambda = [];
for(let i=400;i<1200;i=i+1)
{
    lambda.push(i);
}
let abs = [];
for(let j=0;j<lambda.length;j++)
{
    abs.push(calcAbs(lambda[j]));
}

// Plot
const data = {
    labels: lambda,
    datasets: [{
        label: 'Absorption',
        backgroundColor: 'rgb(255, 99, 132)',
        borderColor: 'rgb(255, 99, 132)',
        pointRadius: 1,
        data: abs
    }]
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