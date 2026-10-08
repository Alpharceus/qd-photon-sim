// Empty-state plate for workspaces that land in a later build step.
import { h } from "./dom.js";
import { icon } from "./icons.js";

export function comingSoon(el, { iconName, title, body, plans = [], action }) {
  el.replaceChildren(
    h("section.ws-soon", { "aria-labelledby": "soon-title" },
      h("div.soon-plate",
        h("div.soon-ico", icon(iconName, { size: 40 })),
        h("h1.soon-title#soon-title", title),
        h("p.soon-body", body),
        plans.length ? h("ul.soon-plans", plans.map((p) => h("li", p))) : null,
        action ? h("a.btn.btn-outline", { href: action.href }, action.label) : null)));
}
