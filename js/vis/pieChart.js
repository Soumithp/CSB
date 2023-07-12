import { color } from "./constant.js";

export default function createPieChart(containerSelector, data) {
  //data process
  data.deaths_age_sex.forEach((d) => {
    d.x = d.x * 1;
    d.y = d.y * 1;
    d.age = d.age * 1;
    d.gender = d.gender * 1;
  });

  //size
  const margin = { top:-50, right: 40, bottom: 40, left: 40 };
  let width, height, innerWidth, innerHeight;

  //create element
  const container = d3.select(containerSelector);
  const svg = container.append("svg");

  //
  resized();
  window.addEventListener("resize", resized);

  return {};

  function resized() {
    width = container.node().clientWidth;
    innerWidth = width - margin.left - margin.right;

    height = container.node().clientHeight;
    innerHeight = height - margin.top - margin.bottom;

    update();
  }

  function update() {
    //svg size
    svg
      .attr("width", width)
      .attr("viewBox", [
        -margin.left - innerWidth / 2,
        -margin.top - innerHeight / 2,
        width,
        height,
      ])
      .attr("overflow", "visible");

    //data
    let angle = 0;
    const dataAgg = d3
      .rollups(
        data.deaths_age_sex,
        (v) => v.length,
        (d) => d.gender
      )
      .map((d) => {
        const rst = {
          gender: d[0],
          deaths: d[1],
          deathsPP: d[1] / data.deaths_age_sex.length,
          startAngle: angle,
          endAngle: angle + (2 * Math.PI * d[1]) / data.deaths_age_sex.length,
        };
        angle = rst.endAngle;
        return rst;
      });

    //arc
    const radius = Math.min(innerWidth, innerHeight) / 2;
    const pieGenerator = d3
      .arc()
      .innerRadius(radius - 40)
      .outerRadius(radius)
      .startAngle((dd) => dd.startAngle)
      .endAngle((dd) => dd.endAngle)
      .padAngle((2 * Math.PI * 1) / 360);

    //mark
    const pie = svg.selectAll("g.pie").data(dataAgg, (d) => d.gender);
    const pieEnter = pie.enter().append("g").attr("class", "pie");
    const pieUpdate = pie.merge(pieEnter);
    pie.exit().remove();

    pieEnter.append("path");
    pieEnter.append("text").attr("class", "label");
    pieEnter.append("text").attr("class", "value");

    pieUpdate
      .select("path")
      .attr("stroke", "none")
      .attr("fill", (d) => color(d.gender))
      .attr("d", pieGenerator);

    pieUpdate
      .select("text.label")
      .attr("font-size", 10)
      .attr("y", -5 - 2)
      .attr("transform", (d) => `translate(${pieGenerator.centroid(d)})`)
      .text((d) => `${d.gender ? "Female" : "Male"}`);

    pieUpdate
      .select("text.value")
      .attr("font-size", 10)
      .attr("y", 5 + 2)
      .attr("transform", (d) => `translate(${pieGenerator.centroid(d)})`)
      .text((d) => `${d3.format(".1%")(d.deathsPP)}`);

    //event
    pieUpdate
      .on("mouseover", (_, d) => {
        d3.selectAll(".dot")
          .classed("highlight", (dd) => dd.gender == d.gender)
          .classed("lowlight", (dd) => dd.gender !== d.gender);
      })
      .on("mouseout", () => {
        d3.selectAll(".dot")
          .classed("highlight", false)
          .classed("lowlight", false);
      });
  }
}
