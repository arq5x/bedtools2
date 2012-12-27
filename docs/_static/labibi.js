// Store editor pop-up help state in localStorage
// so it does not re-pop-up itself between page loads.
// Do not even to pretend to support IE gracefully.
(function($) {

    $(document).ready(function() {
        var box = $("#editor-trap");
        var klass = "toggled";
        var storageKey = "toggled";

        function toggle() {
            box.toggleClass(klass);
            // Store the toggle status in local storage as "has value string" or null
            window.localStorage.setItem(storageKey, box.hasClass(klass) ? "toggled" : "not-toggled");
        }

        box.click(toggle);

        // Check the persistent state of the editor pop-up
        // Note that localStorage does not necessarily support boolean values (ugh!)
        // http://stackoverflow.com/questions/3263161/cannot-set-boolean-values-in-localstorage
        var v = window.localStorage.getItem(storageKey);
        if(v == "toggled" || !v) {
          box.addClass(klass);
        }

    });

})(jQuery);

