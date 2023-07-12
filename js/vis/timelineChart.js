export default function createTimelineChart(containerSelector, data) {
  let onDateSelectedHandler = () => {};
  let brushDateRange;

  //data process
  data.deathdays.forEach((d) => {
    d.dateStr = d.date;
    d.date = new Date(d.date);
    d.deaths = d.deaths * 1;
  });

  //size
  const margin = { top: 20, right: 40, bottom: 40, left: 40 };
  let width, height, innerWidth, innerHeight;

  //scale
  const xScale = d3
    .scaleTime()
    .domain(d3.extent(data.deathdays, (d) => d.date));

  const yScale = d3
    .scaleLinear()
    .domain([0, d3.max(data.deathdays, (d) => d.deaths)])
    .nice();

  //create element
  const container = d3.select(containerSelector);
  const svg = container.append("svg");
  svg.append("g").attr("class", "mark");
  svg.append("g").attr("class", "axis x-axis");
  svg.append("g").attr("class", "axis y-axis");

  //
  resized();
  window.addEventListener("resize", resized);

  return { onDateSelected };

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

    //scale
    xScale.range([0, innerWidth]);
    yScale.range([innerHeight, 0]);

    //mark
    svg
      .select(".mark")
      .selectAll("path.timeline")
      .data([data.deathdays])
      .join("path")
      .attr("class", "timeline")
      .attr("d", (d) =>
        d3
          .line()
          .x((dd) => xScale(dd.date))
          .y((dd) => yScale(dd.deaths))(d)
      );

    //x-axis
    svg
      .select(".x-axis")
      .attr("transform", `translate(${0},${innerHeight})`)
      .call(d3.axisBottom(xScale).tickFormat(d3.timeFormat("%b%d")));

    //y-axis
    svg.select(".y-axis").call(d3.axisLeft(yScale));

    //event
    svg
      .on("mousemove", (e) => {
        const [x, y] = d3.pointer(e);

        const isInArea =
          0 <= x && x <= innerWidth && 0 <= y && y <= innerHeight;

        if (isInArea) {
          const date = invert(x);
          onDateSelectedHandler([date, date]);
        } else {
          onDateSelectedHandler(brushDateRange);
        }
      })
      .on("mouseout", () => {
        onDateSelectedHandler(brushDateRange);
      });

    //brush
    const brushed = ({ selection, sourceEvent }) => {
      const isResized = !sourceEvent;
      if (isResized || !selection) return;

      brushDateRange = selection.map((d) => xScale.invert(d));

      onDateSelectedHandler(brushDateRange);
    };

    const brush = d3
      .brushX()
      .extent([
        [0, 0],
        [innerWidth, innerHeight],
      ])
      .on("end", brushed);

    svg.call(brush);

    if (brushDateRange)
      this.svg.call(
        brush.move,
        brushDateRange.map((d) => xScale(d))
      );
  }

  function invert(x) {
    const date = xScale.invert(x);

    const year = date.getFullYear();
    const month = date.getMonth();
    const day = date.getDate();

    const dateDay = new Date(year, month, day);
    const offsetDay =
      (date.getTime() - dateDay.getTime()) / 3600000 < 12 ? 0 : 1;

    return new Date(year, month, day + offsetDay);
  }

  function onDateSelected(handler) {
    onDateSelectedHandler = typeof handler == "function" ? handler : () => {};
  }
}
