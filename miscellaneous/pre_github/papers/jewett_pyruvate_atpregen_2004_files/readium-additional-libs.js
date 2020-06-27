/**
 * findAndReplaceDOMText v 0.4.6
 * @author James Padolsey http://james.padolsey.com
 * @license http://unlicense.org/UNLICENSE
 *
 * Matches the text of a DOM node against a regular expression
 * and replaces each match (or node-separated portions of the match)
 * in the specified element.
 */


!function(e,t){"object"==typeof module&&module.exports?module.exports=t():e.findAndReplaceDOMText=t()}(this,function(){function e(e){return String(e).replace(/([.*+?^=!:${}()|[\]\/\\])/g,"\\$1")}function t(){return n.apply(null,arguments)||r.apply(null,arguments)}function n(e,n,i,o,d){if(n&&!n.nodeType&&arguments.length<=2)return!1;var a="function"==typeof i;a&&(i=function(e){return function(t,n){return e(t.text,n.startIndex)}}(i));var s=r(n,{find:e,wrap:a?null:i,replace:a?i:"$"+(o||"&"),prepMatch:function(e,t){if(!e[0])throw"findAndReplaceDOMText cannot handle zero-length matches";if(o>0){var n=e[o];e.index+=e[0].indexOf(n),e[0]=n}return e.endIndex=e.index+e[0].length,e.startIndex=e.index,e.index=t,e},filterElements:d});return t.revert=function(){return s.revert()},!0}function r(e,t){return new i(e,t)}function i(e,n){var r=n.preset&&t.PRESETS[n.preset];if(n.portionMode=n.portionMode||o,r)for(var i in r)a.call(r,i)&&!a.call(n,i)&&(n[i]=r[i]);this.node=e,this.options=n,this.prepMatch=n.prepMatch||this.prepMatch,this.reverts=[],this.matches=this.search(),this.matches.length&&this.processMatches()}var o="retain",d=document,a={}.hasOwnProperty;return t.NON_PROSE_ELEMENTS={br:1,hr:1,script:1,style:1,img:1,video:1,audio:1,canvas:1,svg:1,map:1,object:1,input:1,textarea:1,select:1,option:1,optgroup:1,button:1},t.NON_CONTIGUOUS_PROSE_ELEMENTS={address:1,article:1,aside:1,blockquote:1,dd:1,div:1,dl:1,fieldset:1,figcaption:1,figure:1,footer:1,form:1,h1:1,h2:1,h3:1,h4:1,h5:1,h6:1,header:1,hgroup:1,hr:1,main:1,nav:1,noscript:1,ol:1,output:1,p:1,pre:1,section:1,ul:1,br:1,li:1,summary:1,dt:1,details:1,rp:1,rt:1,rtc:1,script:1,style:1,img:1,video:1,audio:1,canvas:1,svg:1,map:1,object:1,input:1,textarea:1,select:1,option:1,optgroup:1,button:1,table:1,tbody:1,thead:1,th:1,tr:1,td:1,caption:1,col:1,tfoot:1,colgroup:1},t.NON_INLINE_PROSE=function(e){return a.call(t.NON_CONTIGUOUS_PROSE_ELEMENTS,e.nodeName.toLowerCase())},t.PRESETS={prose:{forceContext:t.NON_INLINE_PROSE,filterElements:function(e){return!a.call(t.NON_PROSE_ELEMENTS,e.nodeName.toLowerCase())}}},t.Finder=i,i.prototype={search:function(){function t(e){for(var d=0,p=e.length;d<p;++d){var l=e[d];if("string"==typeof l){if(o.global)for(;n=o.exec(l);)a.push(s.prepMatch(n,r++,i));else(n=l.match(o))&&a.push(s.prepMatch(n,0,i));i+=l.length}else t(l)}}var n,r=0,i=0,o=this.options.find,d=this.getAggregateText(),a=[],s=this;return o="string"==typeof o?RegExp(e(o),"g"):o,t(d),a},prepMatch:function(e,t,n){if(!e[0])throw new Error("findAndReplaceDOMText cannot handle zero-length matches");return e.endIndex=n+e.index+e[0].length,e.startIndex=n+e.index,e.index=t,e},getAggregateText:function(){function e(r){if(r.nodeType===Node.TEXT_NODE)return[r.data];if(t&&!t(r))return[];var i=[""],o=0;if(r=r.firstChild)do{if(r.nodeType!==Node.TEXT_NODE){var d=e(r);n&&r.nodeType===Node.ELEMENT_NODE&&(!0===n||n(r))?(i[++o]=d,i[++o]=""):("string"==typeof d[0]&&(i[o]+=d.shift()),d.length&&(i[++o]=d,i[++o]=""))}else i[o]+=r.data}while(r=r.nextSibling);return i}var t=this.options.filterElements,n=this.options.forceContext;return e(this.node)},processMatches:function(){var e,t,n,r=this.matches,i=this.node,o=this.options.filterElements,d=[],a=i,s=r.shift(),p=0,l=0,h=0,c=[i];e:for(;;){if(a.nodeType===Node.TEXT_NODE&&(!t&&a.length+p>=s.endIndex?t={node:a,index:h++,text:a.data.substring(s.startIndex-p,s.endIndex-p),indexInMatch:0===p?0:p-s.startIndex,indexInNode:s.startIndex-p,endIndexInNode:s.endIndex-p,isEnd:!0}:e&&d.push({node:a,index:h++,text:a.data,indexInMatch:p-s.startIndex,indexInNode:0}),!e&&a.length+p>s.startIndex&&(e={node:a,index:h++,indexInMatch:0,indexInNode:s.startIndex-p,endIndexInNode:s.endIndex-p,text:a.data.substring(s.startIndex-p,s.endIndex-p)}),p+=a.data.length),n=a.nodeType===Node.ELEMENT_NODE&&o&&!o(a),e&&t){if(a=this.replaceMatch(s,e,d,t),p-=t.node.data.length-t.endIndexInNode,e=null,t=null,d=[],s=r.shift(),h=0,l++,!s)break}else if(!n&&(a.firstChild||a.nextSibling)){a.firstChild?(c.push(a),a=a.firstChild):a=a.nextSibling;continue}for(;;){if(a.nextSibling){a=a.nextSibling;break}if((a=c.pop())===i)break e}}},revert:function(){for(var e=this.reverts.length;e--;)this.reverts[e]();this.reverts=[]},prepareReplacementString:function(e,t,n){var r=this.options.portionMode;return"first"===r&&t.indexInMatch>0?"":(e=e.replace(/\$(\d+|&|`|')/g,function(e,t){var r;switch(t){case"&":r=n[0];break;case"`":r=n.input.substring(0,n.startIndex);break;case"'":r=n.input.substring(n.endIndex);break;default:r=n[+t]||""}return r}),"first"===r?e:t.isEnd?e.substring(t.indexInMatch):e.substring(t.indexInMatch,t.indexInMatch+t.text.length))},getPortionReplacementNode:function(e,t){var n=this.options.replace||"$&",r=this.options.wrap,i=this.options.wrapClass;if(r&&r.nodeType){var o=d.createElement("div");o.innerHTML=r.outerHTML||(new XMLSerializer).serializeToString(r),r=o.firstChild}if("function"==typeof n)return n=n(e,t),n&&n.nodeType?n:d.createTextNode(String(n));var a="string"==typeof r?d.createElement(r):r;return a&&i&&(a.className=i),n=d.createTextNode(this.prepareReplacementString(n,e,t)),n.data&&a?(a.appendChild(n),a):n},replaceMatch:function(e,t,n,r){var i,o,a=t.node,s=r.node;if(a===s){var p=a;t.indexInNode>0&&(i=d.createTextNode(p.data.substring(0,t.indexInNode)),p.parentNode.insertBefore(i,p));var l=this.getPortionReplacementNode(r,e);return p.parentNode.insertBefore(l,p),r.endIndexInNode<p.length&&(o=d.createTextNode(p.data.substring(r.endIndexInNode)),p.parentNode.insertBefore(o,p)),p.parentNode.removeChild(p),this.reverts.push(function(){i===l.previousSibling&&i.parentNode.removeChild(i),o===l.nextSibling&&o.parentNode.removeChild(o),l.parentNode.replaceChild(p,l)}),l}i=d.createTextNode(a.data.substring(0,t.indexInNode)),o=d.createTextNode(s.data.substring(r.endIndexInNode));for(var h=this.getPortionReplacementNode(t,e),c=[],u=0,f=n.length;u<f;++u){var x=n[u],N=this.getPortionReplacementNode(x,e);x.node.parentNode.replaceChild(N,x.node),this.reverts.push(function(e,t){return function(){t.parentNode.replaceChild(e.node,t)}}(x,N)),c.push(N)}var g=this.getPortionReplacementNode(r,e);return a.parentNode.insertBefore(i,a),a.parentNode.insertBefore(h,a),a.parentNode.removeChild(a),s.parentNode.insertBefore(g,s),s.parentNode.insertBefore(o,s),s.parentNode.removeChild(s),this.reverts.push(function(){i.parentNode.removeChild(i),h.parentNode.replaceChild(a,h),o.parentNode.removeChild(o),g.parentNode.replaceChild(s,g)}),g}},t}),define("FindAndReplace",function(){}),define("readium-additional-libs",function(){});
//# sourceMappingURL=readium-additional-libs.js.map