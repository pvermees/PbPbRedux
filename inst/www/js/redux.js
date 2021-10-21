var inputTable = null;
var outputTable = null;

function start() {
    rrpc.initialize();
    inputTable = createDataEntryGrid('input-table', 10, 100);
    shinylight.initialize();
}

function loaddefaults(fn){
    shinylight.call('loadTable', {fn: fn}, null).then(function(result) {
	result.data = result.data.tab;
        shinylight.setGridResult(inputTable, result);
    }).catch(function(reason) {
        shinylight.setElementText('error', reason);
    });
}
