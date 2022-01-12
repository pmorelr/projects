function changeFields(){
    $('.param').hide();
    if ($('#geo').val()=="jet") {
      $('.param'+'.'+$('#model').val()).show();
    } else {
      $('.param'+'.'+$('#geo').val()).show();
    }
};

function openCloseHelp() {
    var x = document.getElementById("help");
    if (x.style.display === "none") {
        x.style.display = "block";
    } else {
        x.style.display = "none";
    }
} 

function bindAll(){
    $('#helpbutton').click(openCloseHelp);
    $('#geo').change(changeFields);
    $('#model').change(changeFields);
};

$(document).ready(function() {
    $('#parametros').append($('<div id="help" style="display: none">').load('/simuladores/erosion_help'));
    bindAll();

    $('#geo').trigger('change');
});
