export default function createBarChart(containerSelector, data) {
  //data process
  data.deaths_age_sex.forEach((d) => {
    d.x = d.x * 1;
    d.y = d.y * 1;
    d.age = d.age * 1;
    d.gender = d.gender * 1;
  });

  //size
  const margin = { top: 20, right: 40, bottom: 40, left: 40 };
  let width, height, innerWidth, innerHeight;

  //scale
  const xScale = d3.scaleBand().domain([0, 1, 2, 3, 4, 5]).padding(0.2);
  const yScale = d3.scaleLinear();

  //create element
  const container = d3.select(containerSelector);
  const svg = container.append("svg");
  svg.append("g").attr("class", "mark");
  svg.append("g").attr("class", "axis x-axis");
  svg.append("g").attr("class", "axis y-axis");

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
      .attr("viewBox", [-margin.left, -margin.top, width, height]);

    //data
    const dataAgg = d3
      .rollups(
        data.deaths_age_sex,
        (v) => v.length,
        (d) => d.age
      )
      .map((d) => ({ age: d[0], deaths: d[1] }));

    //scale
    xScale.range([0, innerWidth]);
    yScale
      .domain([0, d3.max(dataAgg, (d) => d.deaths)])
      .nice()
      .range([innerHeight, 0]);

    //mark
    const total = d3.sum(dataAgg, (d) => d.deaths);
    svg
      .select(".mark")
        .selectAll("rect.bar", "000000")
      .data(dataAgg)
        .join("rect")
        .attr("class", "bar","fill= 000000")
      .attr("x", (d) => xScale(d.age))
      .attr("y", (d) => yScale(d.deaths))
      .attr("width", xScale.bandwidth)
      .attr("height", (d) => innerHeight - yScale(d.deaths))
      .on("mouseover", (_,d) => {
        svg
            .selectAll("text.label")
            
          .data([""])
          .join("text")
          .attr("class", "label")
          .attr("dominant-baseline", "auto")
          .attr("text-anchor", "middle")
            .attr("fill", "red")
          .attr("font-size", 12)
          .attr("x", xScale(d.age) + xScale.bandwidth() / 2)
          .attr("y", yScale(d.deaths)-5)
          .text(`${d.deaths} (${d3.format(".1%")(d.deaths / total)})`);
      })
      .on("mouseout", () => {
        svg.selectAll("text.label").remove();
      });

    //x-axis
    svg
      .select(".x-axis")
      .attr("transform", `translate(${0},${innerHeight})`)
      .call(
        d3
          .axisBottom(xScale)
          .tickFormat(
            (age) => ["0-10", "11-20", "21-40", "41-60", "61-80", "> 80"][age]
          )
      );

    //y-axis
    svg.select(".y-axis").call(d3.axisLeft(yScale));
  }
}
