/* js/loader.js */
document.addEventListener("DOMContentLoaded", function() {
    // 1. Determina la lingua (se il file finisce con ITA o è DidatticaITA, etc)
    const isIta = window.location.pathname.includes("ITA");
    const headerFile = isIta ? 'components/header-it.html' : 'components/header-en.html';
    
    // 2. Carica Header
    fetch(headerFile)
        .then(response => response.text())
        .then(data => {
            document.getElementById('header-placeholder').innerHTML = data;
            setActiveLink();
        });

    // 3. Carica Footer
    fetch('components/footer.html')
        .then(response => response.text())
        .then(data => {
            document.getElementById('footer-placeholder').innerHTML = data;
        });
});

function setActiveLink() {
    // Evidenzia il link attivo nella navbar basandosi sull'URL corrente
    const path = window.location.pathname;
    const links = document.querySelectorAll('.topnav a');
    
    links.forEach(link => {
        // Logica semplice: se l'href del link è contenuto nell'URL della pagina
        if (path.includes(link.getAttribute('href'))) {
            link.classList.add('active');
        }
    });
}