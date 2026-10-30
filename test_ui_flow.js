'use strict';
// Exercise the actual UI setup-completion handler without a browser dependency.
const fs=require('fs'),vm=require('vm'),assert=require('assert');
const source=fs.readFileSync('app_assets/app.js','utf8');
const handler=source.slice(source.indexOf('async function pollSetup()'),source.indexOf('async function startSetup('));
async function check(state,ready,expected){
 const nodes=new Map();const $=id=>{if(!nodes.has(id))nodes.set(id,{classList:{toggle(){}},value:'',textContent:'',hidden:false});return nodes.get(id)};
 let calls=[];const selected={assembly:'original.fasta'};
 const context={$ ,api:async()=>({done:true,state,install:true,result:{ready},log:''}),navigator:{platform:'Win32'},updateEnvironment(){},updateRunAvailability(){},clearInterval(){},error(){},launchAnalysis:async s=>calls.push(s),pendingRun:selected,setupTimer:1,setupActive:true,runtimeReady:false};
 vm.createContext(context);await vm.runInContext(handler+';pollSetup()',context);
 assert.equal(calls.length,expected);if(expected)assert.strictEqual(calls[0],selected);assert.equal(context.pendingRun,null);
}
(async()=>{await check('finished',true,1);await check('failed',false,0);await check('cancelled',false,0);console.log('3 UI setup handoff checks passed');})().catch(e=>{console.error(e);process.exit(1)});
