import createMap from "./vis/map.js";
import createTimelineChart from "./vis/timelineChart.js";
import createPieChart from "./vis/pieChart.js";
import createBarChart from "./vis/barChart.js";

Promise.all([
  d3.json("../data/streets.json"),
  d3.csv("../data/pumps.csv"),
  d3.csv("../data/deathdays.csv"),
  // d3.csv("../data/deaths_age_sex.csv"),
  d3.csv("../data/deaths_age_sex_with_date.csv"),
]).then(ready);

function ready([streets, pumps, deathdays, deaths_age_sex]) {
  //create vis
  const map = createMap("#vis-map", { streets, pumps, deaths_age_sex });
  const timelineChart = createTimelineChart("#vis-timeline-chart", {
    deathdays,
  });
  const pieChart = createPieChart("#vis-gender-chart", { deaths_age_sex });
  const barCart = createBarChart("#vis-age-chart", { deaths_age_sex });

  //event
  timelineChart.onDateSelected((dateRange) => {  
    map.update(dateRange);
  });
}
