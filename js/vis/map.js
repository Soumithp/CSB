import { color } from "./constant.js";

export default function createMap(containerSelector, data) {
  //data process
  data.pumps.forEach((d) => {
    d.x = d.x * 1;
    d.y = d.y * 1;
    d.type = "Pump";
  });

  const landmark = [
    ...data.pumps,
    { x: 11, y: 13, type: "Work House" },
    { x: 13.8798304, y: 12, type: "Brewery" },
  ];

  data.deaths_age_sex.forEach((d) => {
    d.x = d.x * 1;
    d.y = d.y * 1;
    d.age = d.age * 1;
    d.gender = d.gender * 1;
    d.date = new Date(d.date);
  });

  //size
  const margin = { top: 0, right: 0, bottom: 0, left: 0 };
  const width = 600;
  const height = 600;
  const innerWidth = width - margin.left - margin.right;
  const innerHeight = height - margin.top - margin.bottom;

  let dateRange;

  //scale
  const xScale = d3
    .scaleLinear()
    .domain(
      d3.extent(data.streets.map((line) => line.map((pt) => pt.x)).flat())
    )
    .range([0, innerWidth]);

  const yScale = d3
    .scaleLinear()
    .domain(
      d3.extent(data.streets.map((line) => line.map((pt) => pt.y)).flat())
    )
    .range([innerHeight, 0]);

  const typeColor = d3
    .scaleOrdinal()
    .domain(["Pump", "Work House", "Brewery"])
    .range(["black", "green", "hotpink"]);

  //create element
  const container = d3.select(containerSelector);

  container.append("div").attr("class", "legend").call(legend);

  const svg = container
    .append("svg")
    .attr("width", "100%")
    .attr("viewBox", [-margin.left, -margin.top, width, height])
    .attr("transform", `rotate(1)`);

  update();

  return { update };

  function update(_dateRange) {
    if (_dateRange) dateRange = _dateRange;

    //street
    svg
      .selectAll("path.street")
      .data(data.streets)
      .join("path")
      .attr("class", "street")
      .attr("d", (d) =>
        d3
          .line()
          .x((dd) => xScale(dd.x))
          .y((dd) => yScale(dd.y))(d)
      );

    //street name
      const streetNames = [{ name: "BREWER STREET", x: 12, y: 6, rotate: -35 },
          { name: "OXFORD STREET", x: 11, y: 16, rotate: -25 },
          { name: "RECENT STREET", x: 8, y: 11, rotate: 65 },
          { name: "DEAN STREET", x: 17, y: 15, rotate: 65 },
          { name: "SILVER STREET", x: 11, y: 9, rotate: -30 },
          { name: "golden square", x: 11, y: 7.5, rotate: -28 }      ];


    svg
      .selectAll("text.street-name")
      .data(streetNames)
      .join("text")
      .attr("class", "street-name")
      .attr(
        "transform",
        (d) => `translate(${xScale(d.x)},${yScale(d.y)}) rotate(${d.rotate})`
      )
      .text((d) => d.name);

    //pumps
    const landmarkGroup = svg.selectAll("g.landmark").data(landmark);
    const landmarkGroupEnter = landmarkGroup
      .enter()
      .append("g")
      .attr("class", "landmark");
    const landmarkGroupUpdate = landmarkGroup.merge(landmarkGroupEnter);
    landmarkGroup.exit().remove();

    landmarkGroupEnter
      .append("circle")
      .attr("fill", (d) => typeColor(d.type))
      .attr("r", 3);

    landmarkGroupEnter
      .append("text")
      .attr("x", 4)
      .text((d) => d.type);

    landmarkGroupUpdate.attr(
      "transform",
      (d) => `translate(${xScale(d.x)},${yScale(d.y)})`
    );

    //dot
    const dataSelected = data.deaths_age_sex.filter((d) =>
      dateRange == undefined
        ? true
        : dateRange[0].getTime() <= d.date.getTime() &&
          d.date.getTime() <= dateRange[1].getTime()
    );

    const dot = svg.selectAll("g.dot").data(dataSelected);
    const dotEnter = dot.enter().append("g").attr("class", "dot");
    const dotUpdate = dot.merge(dotEnter);
    dot.exit().remove();

    dotEnter.append("circle");
    dotEnter.append("text");

    dotUpdate.attr(
      "transform",
      (d) => `translate(${xScale(d.x)},${yScale(d.y)})`
    );

    dotUpdate
      .select("circle")
      .attr("stroke", (d) => color(d.gender))
      .attr("fill", (d) => color(d.gender))
      .attr("r", 4);

    dotUpdate
      .select("text")
      .attr("visibility", "hidden")
      .attr("x", 4 + 2)
      .text(
        (d) =>
          `Age ${["0-10", "11-20", "21-40", "41-60", "61-80", "> 80"][d.age]}`
      );

    //event
    dotUpdate
      .on("mouseover", function () {
        dotUpdate.classed("highlight", false).classed("lowlight", true);

        d3.select(this).select("text").attr("visibility", "visible");

        d3.select(this)
          .classed("highlight", true)
          .classed("lowlight", false)
          .raise();
      })
      .on("mouseout", function () {
        d3.select(this).select("text").attr("visibility", "hidden");

        dotUpdate.classed("highlight", false).classed("lowlight", false);
      });

    //zoom
    function handleZoom(e) {
      xScale.range([e.transform.applyX(0), e.transform.applyX(innerWidth)]);
      yScale.range([e.transform.applyY(innerHeight), e.transform.applyY(0)]);
      update();
    }

    const zoom = d3
      .zoom()
      .scaleExtent([1, 50])
      .translateExtent([
        [0, 0],
        [width, height],
      ])
      .on("zoom", handleZoom);
    svg.call(zoom);
  }

  function legend(parent) {
    //
    const item = parent.selectAll("span.legend-item").data(typeColor.domain());
    const itemEnter = item.enter().append("span").attr("class", "legend-item");

    itemEnter
      .append("span")
      .attr("class", "legend-symbol")
      .style("background-color", (d) => typeColor(d));
    itemEnter
      .append("span")
      .attr("class", "legend-label")
      .text((d) => d);
  }
}
