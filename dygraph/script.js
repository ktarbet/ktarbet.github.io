


var g = new Dygraph(document.getElementById("div1"),
    "consumer-debt.csv", {
    legend: 'always',
    title: 'Consumer debt',
    showRoller: true,
    rollPeriod: 14,
    customBars: true,
    ylabel: 'Temperature (F)',
  });

